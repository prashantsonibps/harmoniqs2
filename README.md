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

Our optimizer jointly searches:

- maximum Rabi frequency;
- initial detuning;
- final detuning;
- nonlinear sweep exponent;
- atom spacing.

The optimized protocol uses the full device-valid `6 µs` duration. The
amplitude rises smoothly, remains constant through the detuning sweep, and
falls smoothly to zero.

Optimized `C5` parameters:

```text
maximum Ω       = 5.864953 rad/µs
initial δ       = -12.790586 rad/µs
final δ         = 11.736663 rad/µs
sweep exponent  = 1.017324
duration        = 6000 ns
side length     = 5.261884 µm
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
- optimized dense `P_MIS`: `0.983185`;
- Pulser/QuTiP `P_MIS`: `0.983631`;
- modulation-aware `P_MIS`: `0.984146`.

For `K1,3`:

- reference `P_MIS`: `0.907995`;
- optimized Pulser/QuTiP `P_MIS`: `0.992213`;
- modulation-aware `P_MIS`: `0.993861`.

`C5` is selected for submission because its improvement over the corresponding
reference is larger.

## Repository structure

```text
src/harmoniqs/
  common.py          Device limits, scoring, counts, JSON exports
  challenge01.py     Bell-state model, optimizer, and Pulser validation
  challenge02.py     Graph geometry, MIS model, optimizer, and validation

scripts/
  run_cloud.py       Safe Pasqal Cloud validation and QPU submission
  build_submission.py

tests/
  test_common.py
  test_challenges.py

results/
  challenge01/       Bell pulse sequences, controls, and scores
  challenge02/       MIS sequences, parameters, and scores
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

- Hardware job IDs are not yet recorded.
- The simulation is ideal except for optional channel modulation.
- Large-register exact state-vector simulation is exponentially expensive.
- Challenge 3 and its benchmark-instance pipeline will be added separately.
- Hardware results should be compared with the saved modulation-aware
  predictions before final submission.

## Challenge source

This is an independent implementation created from the public specifications
for [A Real Quantum Hackathon](https://github.com/harmoniqs/a-real-quantum-hackathon).
The challenge concept and event belong to their respective organizers; this
repository contains our control strategies, optimization code, simulations,
tests, and generated results.
