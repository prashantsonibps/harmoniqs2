from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from pulser import AnalogDevice, Pulse, Sequence
from pulser_simulation import QutipEmulator
from scipy.linalg import expm
from scipy.optimize import differential_evolution, minimize

from .common import (
    C6_OVER_HBAR,
    OMEGA_MAX,
    deterministic_counts,
    probabilities_from_state,
    register_with_layout,
    save_json,
    save_sequence,
    validate_sequence,
)

SPACINGS_UM = (5.0, 6.5)
REFERENCE_OMEGA = 2 * np.pi
REFERENCE_DURATION_NS = 352
SEGMENTS = 6
SEGMENT_DURATION_NS = 64
ROBUST_SEGMENT_DURATION_NS = 128
ROBUST_CONTROLS = {
    5.0: np.array(
        [
            [8.52664328, 9.44053212],
            [4.13524880, -3.67173847],
            [6.31903239, -6.85851230],
        ]
    ),
    6.5: np.array(
        [
            [11.79944435, -1.56887568],
            [1.93754212, -4.69787030],
            [4.69154936, 4.97833563],
        ]
    ),
}


def symmetric_hamiltonian(
    omega: float, detuning: float, spacing_um: float
) -> np.ndarray:
    coupling = omega / np.sqrt(2)
    interaction = C6_OVER_HBAR / spacing_um**6
    return np.array(
        [
            [0.0, coupling, 0.0],
            [coupling, -detuning, coupling],
            [0.0, coupling, interaction - 2 * detuning],
        ],
        dtype=float,
    )


def evolve_controls(
    controls: np.ndarray,
    spacing_um: float,
    segment_duration_ns: int = SEGMENT_DURATION_NS,
) -> np.ndarray:
    state = np.array([1.0, 0.0, 0.0], dtype=complex)
    dt_us = segment_duration_ns / 1000
    for omega, detuning in np.asarray(controls).reshape(-1, 2):
        hamiltonian = symmetric_hamiltonian(omega, detuning, spacing_um)
        state = expm(-1j * hamiltonian * dt_us) @ state
    return state


def bell_fidelity_symmetric(state: np.ndarray) -> float:
    return float(abs(np.asarray(state).reshape(-1)[1]) ** 2)


def make_sequence(
    spacing_um: float,
    controls: np.ndarray,
    segment_duration_ns: int = SEGMENT_DURATION_NS,
) -> Sequence:
    register = register_with_layout(
        [(0.0, 0.0), (spacing_um, 0.0)],
        ["q0", "q1"],
    )
    sequence = Sequence(register, AnalogDevice)
    sequence.declare_channel("rydberg", "rydberg_global")
    for omega, detuning in np.asarray(controls).reshape(-1, 2):
        sequence.add(
            Pulse.ConstantPulse(
                segment_duration_ns,
                float(omega),
                float(detuning),
                0.0,
            ),
            "rydberg",
        )
    validate_sequence(sequence)
    return sequence


def reference_sequence(spacing_um: float) -> Sequence:
    return make_sequence(
        spacing_um,
        np.array([[REFERENCE_OMEGA, 0.0]]),
        REFERENCE_DURATION_NS,
    )


def simulate_sequence(
    sequence: Sequence, *, with_modulation: bool = False
) -> tuple[float, dict[str, float]]:
    emulator = QutipEmulator.from_sequence(
        sequence,
        evaluation_times=0.01,
        with_modulation=with_modulation,
    )
    final = emulator.run().get_final_state().full().reshape(-1)
    fidelity = float(abs(final[1] + final[2]) ** 2 / 2)
    return fidelity, probabilities_from_state(final, pulser_basis=True)


def optimize_controls(
    spacing_um: float, *, quick: bool = False, seed: int = 29
) -> tuple[np.ndarray, float]:
    initial = np.column_stack(
        (
            np.full(SEGMENTS, np.pi / (np.sqrt(2) * (SEGMENTS * 0.064))),
            np.zeros(SEGMENTS),
        )
    )
    bounds = [
        bound
        for _ in range(SEGMENTS)
        for bound in ((0.0, OMEGA_MAX), (-25.0, 25.0))
    ]

    def objective(flat: np.ndarray) -> float:
        return 1.0 - bell_fidelity_symmetric(
            evolve_controls(flat.reshape(-1, 2), spacing_um)
        )

    result = differential_evolution(
        objective,
        bounds,
        x0=initial.reshape(-1),
        seed=seed,
        maxiter=18 if quick else 65,
        popsize=6 if quick else 8,
        polish=False,
        workers=1,
        updating="immediate",
    )
    polished = minimize(
        objective,
        result.x,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 500, "ftol": 1e-14},
    )
    controls = polished.x.reshape(-1, 2)
    return controls, 1.0 - float(polished.fun)


def run(output_dir: str | Path = "results/challenge01", quick: bool = False) -> dict:
    output = Path(output_dir)
    report: dict[str, object] = {
        "challenge": 1,
        "target": "(|gr> + |rg>) / sqrt(2)",
        "spacings": {},
    }
    for index, spacing in enumerate(SPACINGS_UM):
        baseline = reference_sequence(spacing)
        baseline_fidelity, baseline_probabilities = simulate_sequence(baseline)
        baseline_modulated, _ = simulate_sequence(
            baseline, with_modulation=True
        )

        controls, model_fidelity = optimize_controls(
            spacing, quick=quick, seed=29 + index
        )
        optimized = make_sequence(spacing, controls)
        fidelity, probabilities = simulate_sequence(optimized)
        modulated_fidelity, modulated_probabilities = simulate_sequence(
            optimized, with_modulation=True
        )
        robust = make_sequence(
            spacing,
            ROBUST_CONTROLS[spacing],
            ROBUST_SEGMENT_DURATION_NS,
        )
        robust_fidelity, robust_probabilities = simulate_sequence(robust)
        robust_modulated_fidelity, robust_modulated_probabilities = (
            simulate_sequence(robust, with_modulation=True)
        )

        spacing_key = f"{spacing:.1f}um"
        save_sequence(optimized, output / f"sequence_{spacing_key}.json")
        save_sequence(
            robust, output / f"sequence_{spacing_key}_robust.json"
        )
        save_json(
            {
                "segment_duration_ns": SEGMENT_DURATION_NS,
                "controls": [
                    {"omega_rad_us": float(row[0]), "detuning_rad_us": float(row[1])}
                    for row in controls
                ],
            },
            output / f"waveform_{spacing_key}.json",
        )
        save_json(
            {
                "segment_duration_ns": ROBUST_SEGMENT_DURATION_NS,
                "controls": [
                    {"omega_rad_us": float(row[0]), "detuning_rad_us": float(row[1])}
                    for row in ROBUST_CONTROLS[spacing]
                ],
            },
            output / f"waveform_{spacing_key}_robust.json",
        )
        report["spacings"][spacing_key] = {
            "baseline_fidelity": baseline_fidelity,
            "baseline_modulated_fidelity": baseline_modulated,
            "optimized_model_fidelity": model_fidelity,
            "optimized_fidelity": fidelity,
            "optimized_modulated_fidelity": modulated_fidelity,
            "strictly_beats_baseline": fidelity > baseline_fidelity,
            "robust_fidelity": robust_fidelity,
            "robust_modulated_fidelity": robust_modulated_fidelity,
            "robust_beats_both_baselines": (
                robust_fidelity > baseline_fidelity
                and robust_modulated_fidelity > baseline_modulated
            ),
            "baseline_counts_500": dict(
                deterministic_counts(baseline_probabilities)
            ),
            "optimized_counts_500": dict(
                deterministic_counts(probabilities)
            ),
            "optimized_modulated_counts_500": dict(
                deterministic_counts(modulated_probabilities)
            ),
            "robust_counts_500": dict(
                deterministic_counts(robust_probabilities)
            ),
            "robust_modulated_counts_500": dict(
                deterministic_counts(robust_modulated_probabilities)
            ),
        }
    save_json(report, output / "scores.json")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Optimize Challenge 01.")
    parser.add_argument("--output", default="results/challenge01")
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    report = run(args.output, args.quick)
    for spacing, values in report["spacings"].items():
        print(
            spacing,
            f"baseline={values['baseline_fidelity']:.8f}",
            f"optimized={values['optimized_fidelity']:.8f}",
            f"modulated={values['optimized_modulated_fidelity']:.8f}",
        )


if __name__ == "__main__":
    main()
