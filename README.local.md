# Local hackathon toolkit

This workspace adds executable Pulser baselines and optimizers to the official
challenge specifications. The event repository itself contains specifications,
not starter code.

## Setup

Python 3.12 is used because the scientific stack is more broadly supported than
the machine's default Python 3.14.

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -e '.[dev]'
```

## Simulate and optimize

```bash
XDG_CACHE_HOME="$PWD/.cache" MPLCONFIGDIR="$PWD/.mplconfig" \
  .venv/bin/python -m harmoniqs.challenge01
XDG_CACHE_HOME="$PWD/.cache" MPLCONFIGDIR="$PWD/.mplconfig" \
  .venv/bin/python -m harmoniqs.challenge02
.venv/bin/python scripts/build_submission.py
```

Use `--quick` on either optimizer for a smoke run. Final sequence JSON,
waveform parameters, deterministic 500-shot population summaries, and scores
are written below `results/`.

## Pasqal Cloud

Keep credentials in the shell environment; never put them in source files:

```bash
export PASQAL_USERNAME='your-email'
export PASQAL_PASSWORD='your-password'
export PASQAL_PROJECT_ID='your-project-id'
```

Validate an exported sequence without consuming shots:

```bash
.venv/bin/python scripts/run_cloud.py \
  results/challenge02/sequence_c5.json
```

Submit only after checking the live device and team run budget:

```bash
.venv/bin/python scripts/run_cloud.py \
  results/challenge02/sequence_c5.json --submit
```

The script asks for a literal `YES` before spending the default 500 shots and
prints the batch/job IDs needed in the Discord submission.

## Notes

- The code reads limits from `pulser.AnalogDevice` at runtime.
- Challenge 1 reports both exact Bell fidelity and modulation-aware previews.
- Challenge 2 scores the same 500-shot budget and exports both `K1,3` and `C5`.
- Challenge 3 is not included: the official repository does not ship the
  paper's benchmark instances, and exact 80-atom state-vector simulation is not
  feasible on a laptop.
