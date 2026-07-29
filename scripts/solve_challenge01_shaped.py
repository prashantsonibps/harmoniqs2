#!/usr/bin/env python3
"""
Challenge 01 — Bell state under blockade
Two Rydberg atoms at spacing r, global Ω(t)/δ(t) drive.
Target: |Ψ⁺⟩ = (|gr⟩+|rg⟩)/√2  from  |gg⟩

Strategy:
  - Direct state evolution in the 4D Hilbert space with substep integration
  - Cost = (1 − F_Ψ⁺) + λ · P_rr(T) — explicitly penalise residual |rr⟩
  - Blue detuning (δ < 0) widens the |s⟩–|rr⟩ gap → suppresses double excitation
  - Active de-excitation: shape δ(t) to sweep any leaked |rr⟩ population
    back into {|gr⟩,|rg⟩} via the Ω-coupled |s⟩ channel
"""
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline
from scipy.linalg import expm
import json, time, sys

# ── Device constants ──────────────────────────────────────────────────
C6_HBAR = 865723.0           # rad/µs · µm⁶
Ω_MAX   = 12.566             # rad/µs  (2π × 2 MHz) — device limit
Δ_MAX   = 125.66             # rad/µs  (2π × 20 MHz)
Ω_REF   = 6.283              # rad/µs  (2π × 1 MHz)
T_REF   = 352                # ns

# ── Pauli / number operators (2-atom, 4D) ────────────────────────────
sx = np.array([[0, 1], [1, 0]], dtype=complex)
n_r = np.array([[0, 0], [0, 1]], dtype=complex)   # |r⟩⟨r|

I2 = np.eye(2, dtype=complex)
sx1 = np.kron(sx, I2);  sx2 = np.kron(I2, sx)
n1  = np.kron(n_r, I2); n2  = np.kron(I2, n_r)
nn  = n1 @ n2

# Basis states
ψ_gg = np.array([1, 0, 0, 0], dtype=complex)
ψ_gr = np.array([0, 1, 0, 0], dtype=complex)
ψ_rg = np.array([0, 0, 1, 0], dtype=complex)
ψ_rr = np.array([0, 0, 0, 1], dtype=complex)
ψ_Ψ  = (ψ_gr + ψ_rg) / np.sqrt(2)


def evolve(Ω, Δ, spacing, substeps=10):
    """
    Evolve |gg⟩ under time-dependent Ω(t), Δ(t).
    Returns final state vector (4D complex).
    """
    V = C6_HBAR / spacing ** 6
    N = len(Ω)
    dt = 1.0 / substeps        # ns per substep  (main step = 1 ns)
    dt_us = dt * 1e-3          # µs per substep

    ψ = ψ_gg.copy()
    for k in range(N):
        H = (Ω[k] / 2) * (sx1 + sx2) - Δ[k] * (n1 + n2) + V * nn
        U = expm(-1j * H * dt_us)
        for _ in range(substeps):
            ψ = U @ ψ
    return ψ


def bell_fidelity(ψ):
    """|⟨Ψ⁺|ψ⟩|²"""
    return abs(np.vdot(ψ_Ψ, ψ)) ** 2


def populations(ψ):
    """Populations [Pgg, Pgr, Prg, Prr]"""
    return np.abs(ψ) ** 2


# ── Parameterisation ──────────────────────────────────────────────────
def decode_params(params, T_ns, n_knots):
    """
    Unbounded params → Ω(t), Δ(t) at 1 ns resolution.
    Uses natural cubic spline through n_knots knots.
    """
    n = n_knots
    # Sigmoid: unbounded → bounded
    Ωk = Ω_MAX * (1.0 / (1.0 + np.exp(-params[:n])))
    Δk = Δ_MAX * (2.0 / (1.0 + np.exp(-params[n:2*n])) - 1.0)
    Ωk = np.clip(Ωk, 0.005, Ω_MAX - 0.005)

    knot_t = np.linspace(0, T_ns, n)
    t_grid = np.arange(T_ns)

    Ω = CubicSpline(knot_t, Ωk, bc_type='natural')(t_grid)
    Δ = CubicSpline(knot_t, Δk, bc_type='natural')(t_grid)

    Ω = np.clip(Ω, 0.0, Ω_MAX)
    Δ = np.clip(Δ, -Δ_MAX, Δ_MAX)

    # Edge apodisation: force smooth ramp-up/down
    # (the spline may bulge at the edges; taper the first/last 10 ns)
    taper = np.ones(T_ns)
    taper[:10] = np.sin(np.linspace(0, np.pi/2, 10))
    taper[-10:] = np.sin(np.linspace(np.pi/2, 0, 10))
    Ω *= taper

    return Ω, Δ, (knot_t, Ωk, Δk)


# ── Cost function ─────────────────────────────────────────────────────
def cost(params, T_ns, n_knots, spacing, λ=0.3, substeps=10):
    """Cost = (1 − F) + λ · P_rr(T)."""
    Ω, Δ, _ = decode_params(params, T_ns, n_knots)
    ψ = evolve(Ω, Δ, spacing, substeps)
    F = bell_fidelity(ψ)
    pops = populations(ψ)
    return (1.0 - F) + λ * pops[3], F, pops[3]


# ── Reference ─────────────────────────────────────────────────────────
def reference(spacing):
    """Reference square-pulse result."""
    Ω = np.full(T_REF, Ω_REF)
    Δ = np.zeros(T_REF)
    ψ = evolve(Ω, Δ, spacing, 10)
    F = bell_fidelity(ψ)
    pops = populations(ψ)
    V_ratio = C6_HBAR / spacing ** 6 / Ω_REF
    return F, pops, V_ratio


# ── Optimizer ─────────────────────────────────────────────────────────
def optimize(spacing, T_ns=1200, n_knots=10, n_restarts=6,
             λ=0.3, max_iter=500):
    """L-BFGS-B with multiple restarts."""
    best = dict(x=None, cost=1e9, F=0, Prr=1)

    def objective(x):
        c, _, _ = cost(x, T_ns, n_knots, spacing, λ=λ)
        return c

    for attempt in range(n_restarts):
        x0 = np.zeros(2 * n_knots)

        # Ω: initial guess ≈ Ω_REF with smooth shape
        p = min(Ω_REF / Ω_MAX, 0.95)
        x0[:n_knots] = np.log(p / (1.0 - p + 1e-10))

        # Δ: start near 0 with slight negative bias (blue detuning)
        # to suppress |rr⟩ during transfer
        for i in range(n_knots):
            frac = i / (n_knots - 1)                     # 0 → 1
            x0[n_knots + i] = -0.3 * np.sin(np.pi * frac)    # U-shaped dip

        if attempt > 0:
            x0 += np.random.randn(2 * n_knots) * 0.5

        _, F0, Prr0 = cost(x0, T_ns, n_knots, spacing, λ=0.0)
        sys.stdout.write(f"  [{attempt+1}/{n_restarts}]  init F={F0:.5f}  Prr={Prr0:.4f}")
        sys.stdout.flush()

        res = minimize(objective, x0, method='L-BFGS-B',
                       options=dict(maxiter=max_iter, ftol=1e-14,
                                    gtol=1e-10, disp=0))

        _, Fb, Prrb = cost(res.x, T_ns, n_knots, spacing, λ=0.0)
        print(f"  →  F={Fb:.10f}  Prr={Prrb:.3e}  (iters={res.nit})")

        if res.fun < best['cost']:
            best = dict(x=res.x.copy(), cost=res.fun, F=Fb, Prr=Prrb)

    return best


# ═════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    np.random.seed(42)

    for spacing, label in [(5.0, 'r1_5um'), (6.5, 'r2_6.5um')]:
        print(f"\n{'='*60}")
        print(f"  r = {spacing} µm  ({label})")
        print(f"{'='*60}")

        t0 = time.time()

        # Reference
        F_ref, pops_ref, V_Ω = reference(spacing)
        print(f"\n  Reference (Ω=Ω_ref, T={T_REF}ns, Δ=0):")
        print(f"    V/Ω = {V_Ω:.1f}")
        print(f"    F = {F_ref:.8f}")
        print(f"    Pops: gg={pops_ref[0]:.4f}  gr={pops_ref[1]:.4f}  "
              f"rg={pops_ref[2]:.4f}  rr={pops_ref[3]:.4f}")

        # Configure per spacing
        if spacing <= 5.0:
            cfg = dict(T_ns=800, n_knots=8, n_restarts=5, λ=0.2)
        else:
            cfg = dict(T_ns=1500, n_knots=10, n_restarts=6, λ=0.4)

        print(f"\n  Optimising: T={cfg['T_ns']}ns, "
              f"knots={cfg['n_knots']}, λ_leak={cfg['λ']}")

        best = optimize(spacing, **cfg)

        Ω_opt, Δ_opt, knots = decode_params(best['x'], cfg['T_ns'],
                                              cfg['n_knots'])
        ψ_final = evolve(Ω_opt, Δ_opt, spacing, substeps=20)
        F_final = bell_fidelity(ψ_final)
        pops_final = populations(ψ_final)

        wall = time.time() - t0

        print(f"\n  === Optimised ===")
        print(f"  F = {F_final:.12f}")
        print(f"  |rr⟩ = {pops_final[3]:.3e}")
        print(f"  Pops: gg={pops_final[0]:.3e}  gr={pops_final[1]:.8f}  "
              f"rg={pops_final[2]:.8f}  rr={pops_final[3]:.3e}")
        print(f"  Ω: [{Ω_opt.min():.3f}, {Ω_opt.max():.3f}] rad/µs")
        print(f"  Δ: [{Δ_opt.min():.3f}, {Δ_opt.max():.3f}] rad/µs")
        print(f"  T = {cfg['T_ns']} ns")
        print(f"  ΔF vs ref = {F_final - F_ref:+.6f}")
        print(f"  Wall: {wall:.1f}s")

        # Save
        data = dict(
            spacing_um=spacing, T_ns=cfg['T_ns'],
            n_knots=cfg['n_knots'],
            fidelity=float(F_final), Prr=float(pops_final[3]),
            reference_fidelity=float(F_ref),
            improvement=float(F_final - F_ref),
            populations={k: float(v) for k, v in
                         zip("gg gr rg rr".split(), pops_final)},
            reference_populations={k: float(v) for k, v in
                                   zip("gg gr rg rr".split(), pops_ref)},
            pulse=dict(
                amplitude_rad_per_us=[float(x) for x in Ω_opt],
                detuning_rad_per_us=[float(x) for x in Δ_opt],
                time_ns=list(range(cfg['T_ns'])),
                Ω_bounds=[float(Ω_opt.min()), float(Ω_opt.max())],
                Δ_bounds=[float(Δ_opt.min()), float(Δ_opt.max())],
            ),
            V_over_Omega=float(V_Ω), wall_seconds=wall,
        )

        path = f"results_{label}.json"
        with open(path, 'w') as f:
            json.dump(data, f, indent=2, default=str)
        print(f"  → {path}")
