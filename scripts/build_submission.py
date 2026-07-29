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
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
