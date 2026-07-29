# Harmoniqs hackathon submission

Highest challenge attempted: **Challenge 01 — Shaped-pulse optimization (2nd pass)**

## Challenge 01 result — shaped pulse (best)

Method: B-spline shaped pulse with direct Hamiltonian integration and explicit leakage penalty.

| Spacing | Duration | $F$ (Pulser @ 4 ns) | $1-F$ | Baseline $F$ | Improvement |
|--------:|---------:|---------------------:|------:|-------------:|-----------:|
| 5.0 µm  | 800 ns   | 0.99999967           | $3.3\times 10^{-7}$ | 0.99256 | $\times 22\,000$ |
| 6.5 µm  | 1500 ns  | 0.99999994           | $6.4\times 10^{-8}$ | 0.75027 | $\times 3.9\times 10^{6}$ |

- Cross-validated: QuTiP and Pulser (QutipBackend) agree to $\Delta F < 2\times 10^{-7}$ at 4 ns clock.
- Both pulses within AnalogDevice bounds ($\Omega_\text{max} = 12.57$ rad/µs, $|\delta|_\text{max} = 125.7$ rad/µs).
- Pulser Studio sequences (clean JSON, not double-escaped): `results/challenge01/shaped/sequence_*.json`
- Optimizer reference: `scripts/solve_challenge01_shaped.py`

## Key insight

At weak blockade ($r_2 = 6.5$ µm, $V/\Omega \approx 1.8$), the reference square pulse leaks 22% into $|rr\rangle$. The shaped pulse starts $\delta(t)$ far blue-detuned ($-46$ rad/µs), widening the $|s\rangle$–$|rr\rangle$ gap to $V - \delta \approx 57$ rad/µs to suppress double excitation, then sweeps through resonance mid-pulse to actively drain residual $|rr\rangle$ via the $\Omega$-coupled $|s\rangle$ channel. This reduces $|rr\rangle$ population from 22% to $< 10^{-7}$.

## Challenge 01 supporting result — segment-based optimization

- 5.0 µm: $F = 0.99999824$ (baseline $0.99256$); robust modulated $F = 0.99850$
- 6.5 µm: $F = 0.99996842$ (baseline $0.75000$); robust modulated $F = 0.99271$

## Challenge 02 result
- Graph: `c5`
- Simulated $P_\text{MIS}$: **0.983631**
- Baseline $P_\text{MIS}$: **0.787583**
- Shots: 500
- Sequence/register: `results/challenge02/sequence_c5.json`
- Pulse parameters: `results/challenge02/parameters_c5.json`
- Pasqal Cloud job IDs: PENDING HARDWARE RUN
