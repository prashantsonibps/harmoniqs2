# Challenge 1 — Bell-State Preparation Under Imperfect Blockade

Prepare the symmetric Bell state on two neutral atoms driven by a **global**
Rydberg laser, at two different atom spacings, and beat the reference square
pulse at both.

## Problem

Two atoms start in `|gg>`. The target is

```text
|Ψ+> = (|gr> + |rg>) / √2
```

Only global controls are available — a single `Ω(t)` and a single `δ(t)` act on
both atoms at once — so the Bell state must be produced by the interplay of the
drive and the Rydberg interaction, not by addressing atoms individually.

The two required spacings differ in how much help the blockade gives:

| Spacing | Regime | V/Ω | What limits fidelity |
|---|---|---|---|
| `5.0 µm` | strong blockade | ≈ 8.8 | `|rr>` is already suppressed; small residual leakage |
| `6.5 µm` | weak blockade | ≈ 1.8 | `|rr>` is only partly suppressed and must be actively cancelled |

**Reference protocol:** a resonant square pulse, `Ω = 2π rad/µs`, `δ = 0`,
`352 ns`. It reaches `0.9926` at 5.0 µm but only `0.7500` at 6.5 µm — the weak
blockade case is where the real work is.

## Physics

In the symmetric subspace the two-atom dynamics close on three states,
`|gg>`, `|Ψ+>`, `|rr>`, with Hamiltonian

```text
H/ℏ = [ 0        Ω/√2      0            ]
      [ Ω/√2     -δ        Ω/√2         ]
      [ 0        Ω/√2      V - 2δ       ] ,   V = C₆ / r⁶
```

Because `|Ψ+>` is the only symmetric state coupling `|gg>` to `|rr>`, the whole
control problem is a three-level population-transfer problem — cheap to
simulate, so the optimizer can afford a genuinely global search.

Implemented in [`symmetric_hamiltonian`](../../src/harmoniqs/challenge01.py#L47).

## Approach

Two independent solutions are in the repository. Both are exported and scored.

### 1. Piecewise-constant segments (main pipeline)

[`src/harmoniqs/challenge01.py`](../../src/harmoniqs/challenge01.py)

- The pulse is **6 independently controlled segments of 64 ns** (384 ns total),
  each with its own `(Ω, δ)`, bounded by the device limits.
- Each segment is propagated exactly with a matrix exponential in the 3-state
  basis; the objective is `1 − |<Ψ+|ψ(T)>|²`.
- **Differential evolution** does the global search (seeded: `29 + index`),
  then **bounded L-BFGS-B** polishes to machine precision.
- The winning controls are rebuilt as a real `pulser.Sequence` of
  `ConstantPulse` segments, validated against `AnalogDevice`, and re-simulated
  with `pulser-simulation` / QuTiP — so the reported number comes from the same
  simulator the hardware pipeline uses, not from the reduced model.
- Every sequence is simulated **twice**: ideal, and with `with_modulation=True`
  so that finite channel bandwidth smears the sharp segment edges.

**Robust variant.** The high-fidelity solution uses fast edges that modulation
degrades (`0.99997 → 0.9696` at 6.5 µm). A second family — **3 segments of
128 ns**, hard-coded in `ROBUST_CONTROLS` — trades a little ideal fidelity for
controls the hardware can actually reproduce, and is the one that beats the
reference under *both* ideal and modulated simulation.

### 2. B-spline shaped pulses (cross-validated)

[`scripts/solve_challenge01_shaped.py`](../../scripts/solve_challenge01_shaped.py)

A smooth-waveform alternative: `Ω(t)` and `δ(t)` are B-splines over 8–10 knots,
integrated directly, with a cost term that **actively de-excites `|rr>`** rather
than merely avoiding it. Longer pulses (`800 ns` at 5.0 µm, `1500 ns` at
6.5 µm) buy near-perfect suppression. Each result is cross-validated across four
independent evaluations — QuTiP and Pulser, at 1 ns and 4 ns time steps — and
agreement to ~1e-7 is what makes the numbers trustworthy.

## Results

Ideal / modulation-aware fidelities from
[`results/challenge01/scores.json`](../../results/challenge01/scores.json):

| Spacing | Reference | Optimized (ideal) | Optimized (modulated) | Robust (ideal) | Robust (modulated) |
|---|---|---|---|---|---|
| 5.0 µm | 0.99256399 | 0.99999824 | 0.99452721 | 0.99799308 | **0.99850463** |
| 6.5 µm | 0.75000360 | 0.99996842 | 0.96961698 | 0.99174364 | **0.99270641** |

Shaped-pulse family
([`results/challenge01/shaped/scores.json`](../../results/challenge01/shaped/scores.json)):

| Spacing | T | Knots | Reference | QuTiP @4 ns | Pulser @4 ns |
|---|---|---|---|---|---|
| 5.0 µm | 800 ns | 8 | 0.99261010 | 0.99999967 | 0.99999967 |
| 6.5 µm | 1500 ns | 10 | 0.75027226 | 0.99999984 | 0.99999994 |

The residual `|rr>` population at 6.5 µm falls to `~7e-9` — the leakage that
costs the reference a quarter of its fidelity is essentially eliminated.

**Deterministic 500-shot summary at 6.5 µm** shows the same story in counts:
reference `{10: 188, 01: 187, 11: 112, 00: 13}` → optimized `{10: 250, 01: 250}`.

## Hardware

Both spacings were submitted to `FRESNEL_CAN1` at 500 shots
([`results/hardware_jobs.json`](../../results/hardware_jobs.json)):

- `6.5 µm` — job `ae3bac5d-…`, **completed**;
- `5.0 µm` — job `cadf5a56-…`, **failed: missing calibration** on the device
  side, not a sequence-validity failure.

## Reproducing

```bash
XDG_CACHE_HOME="$PWD/.cache" MPLCONFIGDIR="$PWD/.mplconfig" \
  .venv/bin/python -m harmoniqs.challenge01
```

Add `--quick` for a short smoke run (fewer DE iterations, smaller population),
or `--output DIR` to write elsewhere. The shaped variant:

```bash
.venv/bin/python scripts/solve_challenge01_shaped.py
```

## Outputs

```text
results/challenge01/
  scores.json                    all fidelities + 500-shot count summaries
  sequence_<spacing>.json        max-ideal-fidelity Pulser sequence
  sequence_<spacing>_robust.json modulation-tolerant sequence
  waveform_<spacing>*.json       raw (Ω, δ) per segment
  shaped/
    scores.json                  B-spline results
    cross_validation_results.json  QuTiP/Pulser × 1 ns/4 ns agreement
    sequence_<spacing>.json
    submission.json
```

## Key design decisions

- **Optimize in the reduced 3-state basis, verify in the full simulator.** The
  symmetric model makes global optimization affordable; Pulser/QuTiP is what
  decides whether the answer is real.
- **Report modulated fidelity, not just ideal.** An 0.99997 pulse that collapses
  to 0.97 under finite bandwidth is not the better pulse.
- **Ship both a peak and a robust family** rather than picking one, so the
  hardware run can use the schedule that survives the channel.

## Caveats

- The simulation is closed-system: no spontaneous emission, laser phase noise,
  atom loss, or readout error. Hardware fidelities will be lower.
- `ROBUST_CONTROLS` are hard-coded results of a prior search, not re-optimized
  on each run — rerunning reproduces the main family, not the robust one.
- Deterministic 500-shot counts are a reporting aid, not sampled statistics.
