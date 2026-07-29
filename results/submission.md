# Harmoniqs hackathon submission

Highest challenge attempted: **Challenge 02**

## Challenge 02 result
- Graph: `c5`
- Simulated P_MIS: **0.983631**
- Baseline P_MIS: **0.787583**
- Shots: 500
- Sequence/register: `results/challenge02/sequence_c5.json`
- Pulse parameters: `results/challenge02/parameters_c5.json`
- Pasqal Cloud job IDs: PENDING HARDWARE RUN

## Challenge 01 supporting result
- 5.0um: F=0.99999824 (baseline 0.99256399); robust modulation preview=0.99850463
- 6.5um: F=0.99996842 (baseline 0.75000360); robust modulation preview=0.99270641

## What changed and why
We replaced the linear baseline sweep with a device-valid nonlinear detuning schedule and jointly tuned the drive, endpoints, and register spacing.
The slower 6 µs evolution allocates more time around the small-gap region while the shaped sweep reduces diabatic transitions.
Joint geometry-and-pulse optimization increased the probability mass on maximum independent sets while preserving the intended unit-disk graph.
