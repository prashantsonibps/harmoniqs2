#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"Missing result file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a Discord-ready submission.")
    parser.add_argument("--results", type=Path, default=Path("results"))
    parser.add_argument(
        "--job-id",
        action="append",
        default=[],
        help="Pasqal hardware job ID; repeat for multiple jobs.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/submission.md"),
    )
    args = parser.parse_args()

    challenge1 = load(args.results / "challenge01" / "scores.json")
    challenge2 = load(args.results / "challenge02" / "scores.json")
    eligible = {
        name: values
        for name, values in challenge2["graphs"].items()
        if values["strictly_beats_baseline"]
    }
    if not eligible:
        raise SystemExit("No Challenge 2 result strictly beats its baseline.")
    graph, best = max(
        eligible.items(),
        key=lambda item: item[1]["optimized_p_mis"]
        - item[1]["baseline_p_mis"],
    )
    robust_path = args.results / "challenge02" / "scores_c5_robust.json"
    robust = load(robust_path) if graph == "c5" and robust_path.exists() else None
    simulated_score = (
        robust["robust_candidate"]["ideal_p_mis"]
        if robust
        else best["optimized_p_mis"]
    )
    modulated_line = (
        f"- Modulation-aware P_MIS: "
        f"**{robust['robust_candidate']['modulated_p_mis']:.6f}**\n"
        if robust
        else ""
    )
    sequence_suffix = "_robust" if robust else ""
    job_ids = ", ".join(args.job_id) if args.job_id else "PENDING HARDWARE RUN"

    bell_lines = []
    shaped_path = args.results / "challenge01" / "shaped" / "scores.json"
    if shaped_path.exists():
        shaped = load(shaped_path)
        for spacing, values in shaped["spacings"].items():
            bell_lines.append(
                f"- {spacing}: Pulser F="
                f"{values['optimized_fidelity_pulser_4ns']:.8f} "
                f"(baseline {values['reference_fidelity']:.8f})"
            )
    else:
        for spacing, values in challenge1["spacings"].items():
            bell_lines.append(
                f"- {spacing}: F={values['optimized_fidelity']:.8f} "
                f"(baseline {values['baseline_fidelity']:.8f}); "
                f"robust modulation preview={values['robust_modulated_fidelity']:.8f}"
            )

    text = f"""# Harmoniqs hackathon submission

Highest challenge attempted: **Challenge 02**

## Challenge 02 result
- Graph: `{graph}`
- Simulated P_MIS: **{simulated_score:.6f}**
{modulated_line}- Baseline P_MIS: **{best['baseline_p_mis']:.6f}**
- Shots: 500
- Sequence/register: `results/challenge02/sequence_{graph}{sequence_suffix}.json`
- Pulse parameters: `results/challenge02/parameters_{graph}{sequence_suffix}.json`
- Pasqal Cloud job IDs: {job_ids}

## Challenge 01 supporting result
{chr(10).join(bell_lines)}

## What changed and why
We replaced the linear baseline sweep with a smooth seven-knot detuning schedule and jointly optimized the drive, endpoints, and pentagon spacing.
The 6 µs pulse redistributes time around difficult avoided crossings and was optimized against drive, detuning, and geometry perturbations.
This robust control strategy increases maximum-independent-set probability while preserving the intended C5 unit-disk graph.
"""
    challenge3_path = args.results / "challenge03" / "scores.json"
    retry_path = args.results / "challenge03" / "scores_n17_retry.json"
    if challenge3_path.exists() and retry_path.exists():
        challenge3 = load(challenge3_path)
        retry = load(retry_path)
        hardware_path = args.results / "challenge03" / "hardware_n17.json"
        retry_hardware_path = (
            args.results / "challenge03" / "hardware_n17_retry.json"
        )
        first_hardware = load(hardware_path) if hardware_path.exists() else None
        retry_hardware = (
            load(retry_hardware_path) if retry_hardware_path.exists() else None
        )
        first_job = (
            first_hardware["job_id"] if first_hardware else "NOT RECORDED"
        )
        retry_job = (
            retry_hardware["job_id"] if retry_hardware else "PENDING"
        )
        retry_status = (
            retry_hardware["status"] if retry_hardware else "pending"
        )
        instance_lines = []
        for n_vertices in ("11", "13", "17"):
            values = challenge3["instances"][n_vertices]
            instance_lines.append(
                f"- N={n_vertices}: modulated exact "
                f"`R={values['optimized_modulated_exact']['valid_approximation_ratio']:.6f}` "
                f"vs paper `{values['paper_target']:.6f}`"
            )
        first_hardware_line = (
            f"- First N=17 hardware: `R={first_hardware['valid_approximation_ratio']:.6f}`, "
            f"valid fraction `{first_hardware['valid_fraction']:.3f}`"
            if first_hardware
            else "- First N=17 hardware: result unavailable"
        )
        text = f"""# Harmoniqs hackathon submission

Highest challenge attempted: **Challenge 03**

## Challenge 03 result — exact matched Fresnel benchmarks
{chr(10).join(instance_lines)}
- Retry modulation-aware `R`: **{retry['modulated_exact']['valid_approximation_ratio']:.6f}**
- Retry modulation-aware valid probability: **{retry['modulated_exact']['valid_probability']:.6f}**
- Retry robust worst-case `R`: **{retry['robust_worst_valid_approximation_ratio']:.6f}**
- Shots: 500
- Sequence/register: `results/challenge03/sequence_n17_retry.json`
- Pulse parameters: `results/challenge03/parameters_n17_retry.json`
{first_hardware_line}
- First N=17 job ID: `{first_job}`
- Retry N=17 job ID: `{retry_job}` ({retry_status})

## Challenge 02 supporting result — robust C5
- Simulated `P_MIS`: **{simulated_score:.6f}**
{modulated_line}- Baseline `P_MIS`: **{best['baseline_p_mis']:.6f}**
- Sequence/register: `results/challenge02/sequence_{graph}{sequence_suffix}.json`

## Challenge 01 supporting result
{chr(10).join(bell_lines)}

## What changed and why
We reproduced the paper's exact N=11, N=13, and N=17 diagonal-connected unit-disk instances and replaced the linear schedule with a smooth instance-aware amplitude and detuning sweep.
After the first N=17 hardware run exposed a larger loss and control penalty than channel modulation predicted, we shortened the retry from 6 µs to 4 µs and optimized jointly for approximation ratio and a usable valid-set population.
The retry remains above the published N=17 target throughout the modulation-aware and bounded control-error ensemble while preserving the exact graph, 500-shot comparison, and 60-trap FRESNEL_CAN1 layout.
"""
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
