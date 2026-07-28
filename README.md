# A Real Quantum Hackathon

**Challenge Specifications & Starter-Kit Baselines**

Harmoniqs × Pasqal × Microsoft · July 29, 2026 · Microsoft Garage, NYC

- Event page: [harmoniqs.ai/hackathon](https://harmoniqs.ai/hackathon)
- Questions during the day: your Discord team channel · [discord.gg/VJnXXeJda](https://discord.gg/VJnXXeJda)
- The three challenges live as issues: [Challenge 01](https://github.com/harmoniqs/a-real-quantum-hackathon/issues/1) · [Challenge 02](https://github.com/harmoniqs/a-real-quantum-hackathon/issues/2) · [Challenge 03](https://github.com/harmoniqs/a-real-quantum-hackathon/issues/3)

## The model

All three challenges use the same system: a register of $N$ atoms evolving under the analog-mode Rydberg Hamiltonian, driven by a single global pulse,

$$
\frac{H(t)}{\hbar} \;=\; \frac{\Omega(t)}{2}\sum_i \sigma_x^{(i)} \;-\; \delta(t)\sum_i n_i \;+\; \sum_{i<j}\frac{C_6}{r_{ij}^6}\, n_i n_j ,
$$

where $\Omega(t)$ is the Rabi frequency and $\delta(t)$ the detuning of the global drive (shared by every atom), $n_i = |r\rangle\langle r|_i$ counts the Rydberg excitation of atom $i$, and $r_{ij}$ is fixed by the register layout. Your control variables are the waveforms $\Omega(t)$, $\delta(t)$ and the atom positions.

The **blockade radius** $R_b$ is the spacing at which the interaction equals the drive,

$$
\frac{C_6}{R_b^6} = \hbar\,\Omega \quad\Longrightarrow\quad R_b = \left(\frac{C_6}{\hbar\,\Omega}\right)^{1/6}.
$$

Two atoms closer than $R_b$ cannot both be excited. All three challenges run on that fact.

## Device envelope

Baselines are stated against the analog device model in Pulser (`pulser.AnalogDevice`); the same limits are enforced on Pasqal Cloud. **Verify `Device.specs` at runtime** — the published envelope is authoritative, not this table.

| Quantity | Symbol | Value | Notes |
|---|---|---|---|
| Interaction coefficient | $C_6/\hbar$ | 865,723 rad µs⁻¹ µm⁶ | Rydberg level 60 |
| Max Rabi frequency | $\Omega_{\max}$ | 12.57 rad/µs (2π × 2 MHz) | global channel |
| Max \|detuning\| | $\|\delta\|_{\max}$ | 125.7 rad/µs (2π × 20 MHz) | |
| Max sequence duration | $T_{\max}$ | 6 000 ns | hard cap for all challenges |
| Min atom spacing | $r_{\min}$ | 5 µm | |
| Max atoms / max radius | — | 80 atoms / 38 µm from origin | |
| Waveform clock period | — | 4 ns | durations in multiples of 4 ns |
| Max runs per job | — | 2 000 shots | |

## Workflow & ranking

**Workflow.** Iterate in simulation; validate on Pasqal Cloud within your team's hardware-run budget (announced in Discord).

**Ranking.** Highest stage completed wins; ties within a stage are broken by that stage's own score.

## The challenges

| | System | Task | Score |
|---|---|---|---|
| [Challenge 01](challenges/challenge-01.md) | 2 atoms | Bell-state prep under blockade, at two spacings | fidelity $F$ |
| [Challenge 02](challenges/challenge-02.md) | 4–5 atoms | embed a target graph, sweep to its maximum independent set | $P_{\mathrm{MIS}}$ |
| [Challenge 03](challenges/challenge-03.md) | 10–85+ atoms | beat a published benchmark curve at matched instance size | paper's metric / $\mathcal{R}$ |

## Submission format

Per team, in your Discord team channel by **16:30**:

1. the highest challenge attempted, and for it: your waveforms and register (Pulser sequence, JSON, or code that generates them);
2. simulated score ($F$, $P_{\mathrm{MIS}}$, or $\mathcal{R}$) with the baseline value you beat, same shots;
3. Pasqal Cloud job IDs for your hardware validation runs;
4. three sentences: what you changed relative to the baseline, and why it worked.

## Read ahead

- [arXiv:1808.10816](https://arxiv.org/abs/1808.10816) — *Quantum Optimization for Maximum Independent Set Using Rydberg Atom Arrays* — the mapping; blockade radius as the unit-disk edge.
- [arXiv:2403.11931](https://arxiv.org/abs/2403.11931) — *Graph Algorithms with Neutral Atom Quantum Processors* — graph problems in analog mode.
- [arXiv:2511.22967](https://arxiv.org/abs/2511.22967) — *Benchmarking neutral atom-based quantum processors at scale* — the curve to beat.
- [arXiv:2201.02773](https://arxiv.org/abs/2201.02773) — *A Survey of Quantum Computing for Finance* — background for the bonus.
