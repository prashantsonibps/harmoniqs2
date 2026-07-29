# Neutral-Atom Quantum Pulse Optimization

An end-to-end implementation for designing, simulating, scoring, and exporting
global Rydberg pulse schedules for Pasqal neutral-atom hardware.

This repository currently focuses on two problems:

1. preparing a two-atom Bell state under imperfect Rydberg blockade; and
2. solving Maximum Independent Set (MIS) on small unit-disk graphs.

It contains our simulation and optimization code, device-valid Pulser sequence
exports, reproducible scores, tests, and a guarded Pasqal Cloud submission
workflow. Challenge 3 is intentionally left for a later iteration.

## Why this project exists

Neutral-atom quantum processors encode qubits in two atomic energy levels:

- `|g>` is the ground state;
- `|r>` is a highly excited Rydberg state.

A global laser drives transitions between these states. Two controls determine
the evolution:

- `Ω(t)`, the Rabi frequency, controls the drive strength;
- `δ(t)`, the detuning, changes the energetic preference for Rydberg
  excitations.

Two atoms in `|r>` interact with strength proportional to `1 / r⁶`. At short
distance this interaction shifts the doubly excited state out of resonance.
This is the Rydberg blockade: nearby atoms are unlikely to be excited
simultaneously.

The system is modeled by

```text
H/ℏ = Ω(t)/2 Σᵢ σˣᵢ - δ(t) Σᵢ nᵢ
      + Σᵢ<ⱼ [C₆/(ℏ rᵢⱼ⁶)] nᵢnⱼ .
```

Our task is therefore a control problem. We choose atom positions and pulse
waveforms so that the final quantum state places as much probability as
possible on the desired answer.

## Challenge 1: Bell-state preparation

### Problem

Two atoms start in the product state

```text
|gg>.
```

The target is the symmetric Bell state

```text
|Ψ+> = (|gr> + |rg>) / √2.
```

The pulse must work at two spacings:

- `5.0 µm`, where blockade is strong;
- `6.5 µm`, where blockade leakage makes `|rr>` more likely.

The reference is a resonant `352 ns` square pulse with
`Ω = 2π rad/µs` and `δ = 0`.

### Our solution

The dynamics are symmetric under exchange of the two atoms, so optimization
can be performed in the three-state basis

```text
|gg>, |Ψ+>, |rr>.
```

`src/harmoniqs/challenge01.py` constructs the corresponding Hamiltonian and
propagates the state through six independently controlled `64 ns` segments.
Each segment has bounded Rabi frequency and detuning. Differential evolution
performs the global search, followed by bounded L-BFGS-B polishing.

We then rebuild the optimized controls as a real `pulser.Sequence` and verify
the final state with `pulser-simulation`. A second simulation enables channel
modulation so sharp waveform edges are evaluated more realistically.

Two sequence families are exported:

- `sequence_<spacing>.json` maximizes ideal simulated Bell fidelity;
- `sequence_<spacing>_robust.json` uses smoother controls that remain strong
  after hardware bandwidth modulation.

### Results

At `5.0 µm`:

- reference fidelity: `0.99256399`;
- optimized ideal fidelity: `0.99999824`;
- robust ideal fidelity: `0.99799308`;
- robust modulation-aware fidelity: `0.99850463`.

At `6.5 µm`:

- reference fidelity: `0.75000360`;
- optimized ideal fidelity: `0.99996842`;
- robust ideal fidelity: `0.99174364`;
- robust modulation-aware fidelity: `0.99270641`.

The weak-blockade spacing shows the largest gain because pulse shaping actively
suppresses leakage into `|rr>` instead of relying on blockade alone.

## Challenge 2: Maximum Independent Set

### Problem

An independent set is a set of graph vertices with no edge between any chosen
pair. A Maximum Independent Set contains as many vertices as possible.

For a neutral-atom register:

- each atom represents one graph vertex;
- a Rydberg excitation means the vertex is selected;
- atoms within the blockade radius represent connected vertices;
- blockade penalizes selecting both endpoints of an edge.

We support two target graphs:

- the four-vertex star `K1,3`;
- the five-vertex cycle `C5`.

The main submission uses `C5`. Its maximum independent sets are:

```text
00101
01001
01010
10010
10100
```

Each valid answer selects two non-adjacent vertices.

### Geometry

The five atoms are placed at the vertices of a regular pentagon. Adjacent
atoms are close enough to encode cycle edges, while pentagon diagonals remain
outside the effective edge radius.

The optimized pentagon side length is:

```text
5.261884 µm
```

This preserves the intended `C5` connectivity while changing interaction
strengths enough to improve the evolution.

### Pulse strategy

The reference protocol uses:

- `4 µs` total duration;
- a `252 ns` amplitude rise and fall;
- constant drive during the middle interval;
- a linear detuning sweep from negative to positive.

The fallback optimizer jointly searches:

- maximum Rabi frequency;
- initial detuning;
- final detuning;
- nonlinear sweep exponent;
- atom spacing.

The robust optimizer then replaces the single-exponent sweep with seven smooth
detuning control points and evaluates every candidate under:

- nominal controls;
- ±2% Rabi-frequency errors;
- ±0.3 rad/µs detuning offsets;
- ±0.02 µm atom-position offsets;
- ideal and channel-modulated Pulser simulation.

The optimized protocol uses the full device-valid `6 µs` duration. The
amplitude rises smoothly, remains constant through the detuning sweep, and
falls smoothly to zero.

Robust `C5` parameters:

```text
maximum Ω       = 5.261431 rad/µs
initial δ       = -18.061794 rad/µs
final δ         = 19.025851 rad/µs
duration        = 6000 ns
side length     = 5.304414 µm
```

### Simulation and scoring

Five qubits produce `2⁵ = 32` computational basis states. The dense simulator:

1. builds the complete Rydberg Hamiltonian, including every pairwise `1/r⁶`
   interaction;
2. initializes the register in the all-ground state;
3. evolves the 32-dimensional state through midpoint time steps;
4. converts the final amplitudes into bitstring probabilities;
5. sums the probabilities of the five maximum independent sets.

The best dense candidate is reconstructed as a Pulser waveform and evaluated
again with QuTiP. Scoring always uses the same logical convention:
`1` means Rydberg-excited.

### Results

For `C5`:

- reference `P_MIS`: `0.787583`;
- fallback Pulser/QuTiP `P_MIS`: `0.983631`;
- robust Pulser/QuTiP `P_MIS`: `0.998799`;
- robust modulation-aware `P_MIS`: `0.998747`;
- robust worst-case ensemble `P_MIS`: `0.997979`.

For `K1,3`:

- reference `P_MIS`: `0.907995`;
- optimized Pulser/QuTiP `P_MIS`: `0.992213`;
- modulation-aware `P_MIS`: `0.993861`.

`C5` is selected for submission because its improvement over the corresponding
reference is larger.

## Challenge 3: matched benchmark-scale MIS

Challenge 3 reproduces the exact 11-, 13-, and 17-vertex diagonal-connected unit-disk
grid instances from arXiv:2511.22967 and its public benchmark repository. Both
smaller instances have a classically certified MIS size of four, while N=17 has
MIS size six. The paper's Fresnel
targets are reconstructed from the original 500-shot replicates before any
candidate is compared against them.

`src/harmoniqs/challenge03.py` uses a matrix-free second-order split propagator:
the interaction and detuning phases are diagonal, while the global drive
factorizes into one-qubit rotations. This keeps the complete `2^N` statevector
and every pairwise Rydberg interaction while making exact N=17 simulation
practical. The optimized 6 µs schedule uses smooth amplitude and monotone
detuning control points and is checked with Pulser channel modulation, 500-shot
reporting, time-step convergence, and drive, detuning, and spacing variations.

### Challenge 3 results

The same jointly optimized schedule beats the paper's matched Fresnel curve in
modulation-aware exact simulation at all three reproduced sizes:

- N=11: `R=0.998317` versus the published `0.906616`;
- N=13: `R=0.998583` versus the published `0.908177`;
- N=17: `R=0.979461` versus the published `0.870044`.

The first N=17 hardware run used the 6 µs schedule and obtained `R=0.845238`,
a `0.392` valid-set fraction, and `0.144` MIS probability over 500 shots. That
run exposed a larger hardware penalty than channel modulation alone predicted.
The hardware-calibrated retry therefore shortens the pulse to 4 µs and balances
solution quality against validity; it reaches modulation-aware `R=0.980854`,
valid probability `0.620625`, and robust worst-case `R=0.979088`.

Hardware records:

- first N=17 batch `22954a95-1a6a-44cd-944d-7330f40be804`, job
  `d89cc6b8-4f8a-4c17-81f0-483540effdf4` — completed;
- 4 µs retry batch `640eb36a-1713-4573-9888-98890b7ee080`, job
  `46b2678b-5e6f-48f6-ba79-599ce12637b7` — result pending.

Run the fast reproducible optimization:

```bash
XDG_CACHE_HOME="$PWD/.cache" MPLCONFIGDIR="$PWD/.mplconfig" \
  .venv/bin/python -m harmoniqs.challenge03 --quick
```

All Challenge 3 sequences, 60-trap registers, parameters, exact certificates,
score reports, and provenance are written under `results/challenge03/`.
Hardware submission is intentionally excluded from this command.

## Repository structure

```text
src/harmoniqs/
  common.py          Device limits, scoring, counts, JSON exports
  challenge01.py     Bell-state model, optimizer, and Pulser validation
  challenge02.py     Graph geometry, MIS model, optimizer, and validation
  challenge02_robust.py
                     Smooth ensemble-robust C5 optimization
  challenge03.py     Exact N=11/N=13/N=17 benchmark and robust optimizer

scripts/
  run_cloud.py       Safe Pasqal Cloud validation and QPU submission
  build_submission.py

tests/
  test_common.py
  test_challenges.py

results/
  challenge01/       Bell pulse sequences, controls, and scores
  challenge02/       MIS sequences, parameters, and scores
  challenge03/       Matched benchmark instances, sequences, and scores
  submission.md      Discord-ready submission draft
```

## Local setup

Python 3.12 is recommended.

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -e '.[dev]'
```

Run the tests:

```bash
XDG_CACHE_HOME="$PWD/.cache" MPLCONFIGDIR="$PWD/.mplconfig" \
  .venv/bin/python -m pytest -q
```

Run Challenge 1:

```bash
XDG_CACHE_HOME="$PWD/.cache" MPLCONFIGDIR="$PWD/.mplconfig" \
  .venv/bin/python -m harmoniqs.challenge01
```

Run Challenge 2:

```bash
XDG_CACHE_HOME="$PWD/.cache" MPLCONFIGDIR="$PWD/.mplconfig" \
  .venv/bin/python -m harmoniqs.challenge02
```

Use `--quick` for a shorter smoke optimization.

## Pasqal hardware workflow

The repository never stores cloud credentials. Export them only in your shell:

```bash
export PASQAL_USERNAME='your-email'
export PASQAL_PROJECT_ID='your-project-id'
read -s PASQAL_PASSWORD
export PASQAL_PASSWORD
```

First validate the `C5` sequence without spending hardware shots:

```bash
.venv/bin/python scripts/run_cloud.py \
  results/challenge02/sequence_c5.json
```

If the event specifies a device:

```bash
.venv/bin/python scripts/run_cloud.py \
  results/challenge02/sequence_c5.json \
  --device DEVICE_NAME
```

After confirming the team run budget, submit 500 shots:

```bash
.venv/bin/python scripts/run_cloud.py \
  results/challenge02/sequence_c5.json \
  --device DEVICE_NAME \
  --shots 500 \
  --submit
```

The script requires an explicit `YES` before submission and prints the batch
and job IDs.

The Challenge 3 retry can be regenerated and validated without spending shots:

```bash
.venv/bin/python -m harmoniqs.challenge03 --hardware-retry
.venv/bin/python scripts/run_cloud.py \
  results/challenge03/sequence_n17_retry.json \
  --device FRESNEL_CAN1
```

Regenerate the submission after receiving a job ID:

```bash
.venv/bin/python scripts/build_submission.py \
  --job-id YOUR_JOB_ID
```

## Device-safety checks

Sequences are validated against `pulser.AnalogDevice` at runtime. The current
implementation enforces:

- no more than 80 atoms;
- at least `5 µm` between atoms;
- at most `6000 ns` sequence duration;
- `4 ns` waveform clock alignment;
- bounded Rabi frequency and detuning;
- no more than 2000 requested shots.

The cloud runner fetches the currently available hardware device, switches the
sequence to that device in strict mode, and validates the register again before
submission.

## Reproducibility

- Optimization seeds are fixed.
- Baselines and optimized schedules use the same scoring code.
- Saved scores include deterministic 500-shot population summaries.
- Exported sequences use Pulser's abstract JSON representation.
- Tests cover graph embeddings, device limits, Bell fidelity, MIS scoring, and
  sequence serialization.

The deterministic count summaries are reporting aids, not substitutes for
hardware sampling. Real QPU counts will fluctuate and include preparation,
control, readout, and atom-loss errors.

## Current limitations and next work

- Challenge 3 hardware IDs and the completed first N=17 result are recorded;
  the 4 µs retry result is pending.
- The simulation includes optional channel modulation and bounded control and
  geometry perturbations, but not a complete calibrated SPAM model.
- Large-register exact state-vector simulation is exponentially expensive.
- Challenge 3 exact simulation currently targets N=11, N=13, and N=17; larger
  instances require a tensor-network, blockade-subspace, or QPU backend.
- Hardware results should be compared with the saved modulation-aware
  predictions before final submission.

## Challenge source

This is an independent implementation created from the public specifications
for [A Real Quantum Hackathon](https://github.com/harmoniqs/a-real-quantum-hackathon).
The challenge concept and event belong to their respective organizers; this
repository contains our control strategies, optimization code, simulations,
tests, and generated results.
