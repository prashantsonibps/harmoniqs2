# Challenge 2 — Maximum Independent Set on a Neutral-Atom Register

Solve MIS on small unit-disk graphs by encoding the graph into atom positions
and driving the register with a single global pulse, then measuring.

## Problem

An **independent set** is a set of vertices with no edge between any two of
them; a **maximum** independent set is one of the largest such sets.

The neutral-atom encoding is direct:

| Graph object | Physical counterpart |
|---|---|
| vertex | one atom |
| vertex selected | atom in the Rydberg state `|r>` |
| edge | atoms within the blockade radius |
| independence constraint | blockade penalizes exciting both endpoints |

The figure of merit is `P_MIS` — the total measurement probability of all
bitstrings that are maximum independent sets, with the convention `1 = Rydberg
excited`.

Two graphs are supported:

- **`K1,3`** (four-vertex star) — unique MIS `0111` (the three leaves);
- **`C5`** (five-vertex cycle) — five degenerate MIS solutions,
  `00101`, `01001`, `01010`, `10010`, `10100`.

`C5` is the submitted graph: degenerate ground states and a weaker
reference score make it the harder and more informative instance.

## Geometry

`C5` is embedded as a **regular pentagon** of circumradius
`s / (2 sin(π/5))`. This is what makes the encoding exact: pentagon *sides* fall
inside the blockade radius (cycle edges), while *diagonals* — longer by a factor
of the golden ratio, and `1/r⁶` down by ≈ 18× in interaction — fall outside. The
star is embedded as a center atom with three leaves at 120°.

Side length is itself an optimization variable, not a fixed constant: it tunes
the ratio of edge interaction to drive strength while preserving connectivity.

Implemented in [`register_coordinates`](../../src/harmoniqs/challenge02.py#L45).

## Pulse strategy

**Reference protocol:** `4 µs`, `Ω = 2π rad/µs`, `252 ns` rise and fall, linear
detuning sweep `−4π → +4π`, spacing `5.5 µm`. This is adiabatic state
preparation: start deep in the `δ < 0` regime where `|00000>` is the ground
state, sweep to `δ > 0` where the ground state is a superposition of the maximum
independent sets, and go slowly enough to follow.

Two optimizers improve on it.

### Fallback optimizer — [`challenge02.py`](../../src/harmoniqs/challenge02.py)

Jointly searches five parameters: peak `Ω`, initial `δ`, final `δ`, a **sweep
nonlinearity exponent** `power`, and the atom spacing. Duration is pushed to the
device maximum of `6 µs` — more time is more adiabaticity, for free.

The exponent matters: a linear sweep spends most of its time far from the
avoided crossing where nothing happens. `power ≈ 1.6` slows the sweep exactly
where the gap is smallest.

Search is differential evolution followed by Powell polishing, scored on a dense
32-dimensional propagator (`evolve_dense`) that builds the full Hamiltonian with
**every** pairwise `1/r⁶` term — diagonals included, so the optimizer cannot
cheat by pretending the graph is exactly unit-disk.

### Robust optimizer — [`challenge02_robust.py`](../../src/harmoniqs/challenge02_robust.py)

Two changes:

1. **Free-form sweep.** The single exponent is replaced by **7 PCHIP knots**, so
   `δ(t)` can take any monotone-ish shape the optimizer wants.
2. **Ensemble scoring.** Each candidate is scored not once but across a
   perturbation ensemble, and optimized on its **worst case**:
   - `±2%` Rabi-frequency error,
   - `±0.3 rad/µs` detuning offset,
   - `±0.02 µm` atom-position error.

Optimizing the nominal score alone produces a solution perched on a narrow peak.
Optimizing the worst case produces one that survives real hardware — and here it
also happens to score higher nominally.

Final robust `C5` parameters
([`parameters_c5_robust.json`](../../results/challenge02/parameters_c5_robust.json)):

```text
Ω max        = 5.261431 rad/µs
δ initial    = -18.061794 rad/µs
δ final      =  19.025851 rad/µs
δ knots      = [-18.062, -14.469, -13.301, -1.564, 11.320, 14.878, 19.026]
duration     = 6000 ns
side length  = 5.304414 µm
```

## Simulation and scoring

1. Build the full `2⁵ = 32`-dimensional Rydberg Hamiltonian including all
   pairwise interactions.
2. Initialize in the all-ground state.
3. Evolve through midpoint time steps (72 steps; 44 in `--quick`), diagonalizing
   at each step.
4. Convert final amplitudes to bitstring probabilities.
5. Sum the probabilities of the maximum independent sets.

The best dense candidate is then rebuilt as a real Pulser waveform
(`CompositeWaveform` of ramp / hold / ramp amplitude, interpolated detuning),
validated against `AnalogDevice`, and re-scored with QuTiP — both ideally and
with channel modulation.

## Results

From [`scores.json`](../../results/challenge02/scores.json) and
[`scores_c5_robust.json`](../../results/challenge02/scores_c5_robust.json):

### `C5`

| Protocol | P_MIS |
|---|---|
| reference | 0.787583 |
| fallback (Pulser/QuTiP) | 0.983631 |
| robust (Pulser/QuTiP) | 0.998799 |
| robust, modulation-aware | 0.998747 |
| robust, **worst case over ensemble** | 0.997979 |

The worst-case number is the important one: across every perturbation tested,
the protocol never drops below `0.9980`.

### `K1,3`

| Protocol | P_MIS |
|---|---|
| reference | 0.907995 |
| optimized | 0.992213 |
| modulation-aware | 0.993861 |

`C5` is submitted because its improvement over its own reference is far larger
(`+0.211` vs `+0.084`).

**500-shot summary, `C5` optimized:** the five MIS bitstrings take 98 counts
each (490/500), with the remaining 10 spread over single-excitation states.

## Hardware

Robust `C5` ran on `FRESNEL_CAN1`, 500 shots, job `eba8f06a-…`
([`hardware_c5_robust.json`](../../results/challenge02/hardware_c5_robust.json)):

| Metric | Value |
|---|---|
| hardware `P_MIS` | **0.708** |
| valid independent sets | 432 / 500 (0.864) |
| mean size over valid sets | 1.817 |
| approximation ratio over valid sets | 0.909 |
| simulated prediction | 0.9987 |

The gap between `0.9987` simulated and `0.708` measured is the honest headline
of this challenge: it is decoherence, atom loss, and readout error, none of
which the closed-system simulation contains. What *does* survive is the
structure — the five MIS bitstrings are the five most frequent outcomes
(`01010`: 79, `00101`: 77, `01001`: 72, `10010`: 67, `10100`: 59), and 86% of
shots are at least valid independent sets. The physics encoding works; the
fidelity budget is hardware-limited.

The supporting `K1,3` star sequence was also submitted to `FRESNEL_CAN1` for
500 shots (job `3db5310b-…`). Its current status and complete cloud identifiers
are tracked in [`results/hardware_jobs.json`](../../results/hardware_jobs.json).

## Reproducing

```bash
# fallback optimizer, both graphs
XDG_CACHE_HOME="$PWD/.cache" MPLCONFIGDIR="$PWD/.mplconfig" \
  .venv/bin/python -m harmoniqs.challenge02

# one graph only
.venv/bin/python -m harmoniqs.challenge02 --graph c5

# robust ensemble optimizer (C5)
.venv/bin/python -m harmoniqs.challenge02_robust
```

`--quick` shortens every search; `--output DIR` redirects results.

## Outputs

```text
results/challenge02/
  scores.json                 baseline / optimized / modulated P_MIS per graph
  scores_c5_robust.json       robust candidate + full ensemble scores
  parameters_c5.json          fallback C5 parameters
  parameters_c5_robust.json   robust C5 parameters (7 detuning knots)
  parameters_star.json        K1,3 parameters
  sequence_c5.json            device-valid Pulser sequence (fallback)
  sequence_c5_robust.json     device-valid Pulser sequence (robust, submitted)
  sequence_star.json
  hardware_c5_robust.json     measured 500-shot QPU counts
```

## Key design decisions

- **Optimize geometry alongside the pulse.** Spacing is a control, not a
  constant — the same graph has many valid embeddings with very different gaps.
- **Score on the worst case, not the mean.** A protocol is only as good as its
  behaviour under the errors the machine actually has.
- **Keep every `1/r⁶` term.** Diagonal interactions are small but non-zero;
  dropping them would let the optimizer exploit a graph that doesn't exist.
- **Use the full 6 µs.** Adiabaticity is the cheapest resource available.

## Caveats

- Exact state-vector simulation is `2ⁿ` — this approach does not scale past
  small registers without a different simulator.
- The simulation is closed-system; compare against the hardware numbers above,
  not the ideal ones, when estimating what a new instance will do.
- `FALLBACK_PARAMETERS` in `challenge02_robust.py` are hard-coded from a prior
  run and used as the robust optimizer's starting point.
