# Harmoniqs hackathon submission

Highest challenge attempted: **Challenge 02**

## Challenge 02 result
- Graph: `c5`
- Simulated P_MIS: **0.998799**
- Modulation-aware P_MIS: **0.998747**
- Worst-case ensemble P_MIS: **0.997979**
- Baseline P_MIS: **0.787583**
- Shots: 500
- Sequence/register: `results/challenge02/sequence_c5_robust.json`
- Pulse parameters: `results/challenge02/parameters_c5_robust.json`
- Pasqal Cloud batch ID: `9eae6251-45f9-4c82-a45d-db30eee86c8b`
- Pasqal Cloud job ID: `eba8f06a-d32c-4a58-86af-5e8a0d378ef8`

## Challenge 01 supporting result
- 5.0um: F=0.99999824 (baseline 0.99256399); robust modulation preview=0.99850463
- 6.5um: F=0.99996842 (baseline 0.75000360); robust modulation preview=0.99270641

## What changed and why
We replaced the linear baseline sweep with a smooth seven-knot detuning schedule and jointly optimized the drive, endpoints, and pentagon spacing.
The 6 µs pulse redistributes time around difficult avoided crossings and was optimized against drive, detuning, and geometry perturbations.
This robust control strategy raises modulation-aware P_MIS from the 0.787583 baseline to 0.998747 while preserving the intended C5 unit-disk graph.
