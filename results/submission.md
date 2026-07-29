# Harmoniqs hackathon submission

Highest challenge attempted: **Challenge 02**

## Challenge 02 result — robust C5

- Simulated `P_MIS`: **0.998799**
- Modulation-aware `P_MIS`: **0.998747**
- Worst-case ensemble `P_MIS`: **0.997979**
- Baseline `P_MIS`: **0.787583**
- Shots: 500
- Sequence/register: `results/challenge02/sequence_c5_robust.json`
- Pulse parameters: `results/challenge02/parameters_c5_robust.json`
- Pasqal Cloud batch ID: `9eae6251-45f9-4c82-a45d-db30eee86c8b`
- Pasqal Cloud job ID: `eba8f06a-d32c-4a58-86af-5e8a0d378ef8`
- Hardware `P_MIS`: **0.708** (354/500 shots)
- Hardware valid-independent-set fraction: **0.864** (432/500 shots)
- Hardware counts: `results/challenge02/hardware_c5_robust.json`

## Challenge 02 supporting result — star K1,3

- Simulated `P_MIS`: **0.992213**
- Modulation-aware `P_MIS`: **0.993861**
- Baseline `P_MIS`: **0.907995**
- Sequence/register: `results/challenge02/sequence_star.json`
- Pasqal Cloud batch ID: `3ad69bb5-fa23-4228-8d2e-f0f476c5aaf8`
- Pasqal Cloud job ID: `3db5310b-9a4e-4e48-be1e-951b3b134dca`
- Hardware status: submitted

## Challenge 01 supporting result — shaped pulse

- 5.0 µm: `F = 0.99999967` (baseline `0.99256`)
- 6.5 µm: `F = 0.99999994` (baseline `0.75027`)
- Cross-validation difference: `ΔF < 2 × 10⁻⁷`
- Sequences: `results/challenge01/shaped/sequence_*.json`
- Optimizer: `scripts/solve_challenge01_shaped.py`
- Hardware job records: `results/hardware_jobs.json`

## What changed and why

We replaced the linear baseline sweep with a smooth seven-knot detuning schedule and jointly optimized the drive, endpoints, and pentagon spacing.
The 6 µs pulse redistributes time around difficult avoided crossings and was optimized against drive, detuning, and geometry perturbations.
This robust control strategy raises modulation-aware `P_MIS` from the `0.787583` baseline to `0.998747` while preserving the intended C5 unit-disk graph.
