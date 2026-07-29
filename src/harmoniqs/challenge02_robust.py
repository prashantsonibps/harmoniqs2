from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from pulser import AnalogDevice, Pulse, Sequence
from pulser.waveforms import (
    CompositeWaveform,
    ConstantWaveform,
    InterpolatedWaveform,
    RampWaveform,
)
from scipy.interpolate import PchipInterpolator
from scipy.linalg import eigh
from scipy.optimize import differential_evolution

from .challenge02 import (
    GRAPH_EDGES,
    RISE_NS,
    _operators,
    make_register,
    make_sequence,
    register_coordinates,
    simulate_sequence,
)
from .common import (
    mis_bitstrings,
    probabilities_from_state,
    save_json,
    save_sequence,
    validate_sequence,
)

KNOT_TIMES = np.linspace(0.0, 1.0, 7)
FALLBACK_PARAMETERS = {
    "omega": 5.864953121759875,
    "detuning_initial": -12.790585721504595,
    "detuning_final": 11.736663117195056,
    "power": 1.0173239794790947,
    "duration_ns": 6000,
    "spacing_um": 5.261884497707749,
}


def unpack_parameters(values: np.ndarray) -> dict[str, object]:
    omega, initial, final, spacing = values[:4]
    logits = np.append(values[4:], 0.0)
    weights = np.exp(logits - np.max(logits))
    weights /= weights.sum()
    fractions = np.concatenate(([0.0], np.cumsum(weights)))
    knots = initial + (final - initial) * fractions
    return {
        "omega": float(omega),
        "detuning_initial": float(initial),
        "detuning_final": float(final),
        "detuning_knots": knots.tolist(),
        "duration_ns": 6000,
        "spacing_um": float(spacing),
    }


def make_robust_sequence(parameters: dict[str, object]) -> Sequence:
    duration_ns = int(parameters["duration_ns"])
    hold_ns = duration_ns - 2 * RISE_NS
    omega = float(parameters["omega"])
    initial = float(parameters["detuning_initial"])
    final = float(parameters["detuning_final"])
    knots = np.asarray(parameters["detuning_knots"], dtype=float)

    amplitude = CompositeWaveform(
        RampWaveform(RISE_NS, 0.0, omega),
        ConstantWaveform(hold_ns, omega),
        RampWaveform(RISE_NS, omega, 0.0),
    )
    detuning = CompositeWaveform(
        ConstantWaveform(RISE_NS, initial),
        InterpolatedWaveform(hold_ns, knots, times=KNOT_TIMES),
        ConstantWaveform(RISE_NS, final),
    )
    sequence = Sequence(
        make_register("c5", float(parameters["spacing_um"])),
        AnalogDevice,
    )
    sequence.declare_channel("rydberg", "rydberg_global")
    sequence.add(Pulse(amplitude, detuning, 0.0), "rydberg")
    validate_sequence(sequence)
    return sequence


def evolve_robust(
    parameters: dict[str, object],
    *,
    steps: int = 80,
    omega_scale: float = 1.0,
    detuning_offset: float = 0.0,
    spacing_offset: float = 0.0,
) -> np.ndarray:
    spacing = float(parameters["spacing_um"]) + spacing_offset
    coordinates = register_coordinates("c5", spacing)
    drive, excitations, interactions = _operators(coordinates)
    duration_us = float(parameters["duration_ns"]) / 1000
    rise_us = RISE_NS / 1000
    hold_us = duration_us - 2 * rise_us
    knots = np.asarray(parameters["detuning_knots"], dtype=float)
    interpolator = PchipInterpolator(KNOT_TIMES, knots)
    dt = duration_us / steps
    state = np.zeros(drive.shape[0], dtype=complex)
    state[0] = 1.0

    for step in range(steps):
        time_us = (step + 0.5) * dt
        if time_us < rise_us:
            omega = (
                float(parameters["omega"])
                * omega_scale
                * time_us
                / rise_us
            )
            detuning = knots[0]
        elif time_us <= rise_us + hold_us:
            omega = float(parameters["omega"]) * omega_scale
            x = (time_us - rise_us) / hold_us
            detuning = float(interpolator(x))
        else:
            omega = (
                float(parameters["omega"])
                * omega_scale
                * (duration_us - time_us)
                / rise_us
            )
            detuning = knots[-1]
        hamiltonian = omega * drive + np.diag(
            interactions - (detuning + detuning_offset) * excitations
        )
        eigenvalues, eigenvectors = eigh(hamiltonian)
        state = eigenvectors @ (
            np.exp(-1j * eigenvalues * dt)
            * (eigenvectors.conj().T @ state)
        )
    return state


def score_state(state: np.ndarray) -> float:
    targets = mis_bitstrings(5, GRAPH_EDGES["c5"])
    probabilities = probabilities_from_state(state)
    return float(sum(probabilities[target] for target in targets))


def robust_dense_scores(parameters: dict[str, object]) -> list[float]:
    variations = (
        {},
        {"omega_scale": 0.98},
        {"omega_scale": 1.02},
        {"detuning_offset": -0.3},
        {"detuning_offset": 0.3},
        {"spacing_offset": -0.02},
        {"spacing_offset": 0.02},
    )
    return [
        score_state(evolve_robust(parameters, **variation))
        for variation in variations
    ]


def optimize(seed: int = 113, quick: bool = False) -> tuple[dict[str, object], list[float]]:
    bounds = [
        (4.5, 9.5),
        (-20.0, -8.0),
        (8.0, 20.0),
        (5.02, 6.3),
        *[(-2.5, 2.5) for _ in range(5)],
    ]
    initial = np.array(
        [
            FALLBACK_PARAMETERS["omega"],
            FALLBACK_PARAMETERS["detuning_initial"],
            FALLBACK_PARAMETERS["detuning_final"],
            FALLBACK_PARAMETERS["spacing_um"],
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        dtype=float,
    )

    def objective(values: np.ndarray) -> float:
        scores = robust_dense_scores(unpack_parameters(values))
        return -(0.4 * scores[0] + 0.6 * min(scores))

    result = differential_evolution(
        objective,
        bounds,
        x0=initial,
        seed=seed,
        maxiter=8 if quick else 28,
        popsize=4 if quick else 6,
        polish=True,
        workers=1,
        updating="immediate",
        tol=1e-7,
    )
    parameters = unpack_parameters(result.x)
    return parameters, robust_dense_scores(parameters)


def run(output_dir: str | Path = "results/challenge02", quick: bool = False) -> dict:
    output = Path(output_dir)
    fallback = make_sequence("c5", FALLBACK_PARAMETERS)
    fallback_ideal, _ = simulate_sequence(fallback, GRAPH_EDGES["c5"])
    fallback_modulated, _ = simulate_sequence(
        fallback,
        GRAPH_EDGES["c5"],
        with_modulation=True,
    )

    parameters, dense_scores = optimize(quick=quick)
    candidate = make_robust_sequence(parameters)
    candidate_ideal, _ = simulate_sequence(candidate, GRAPH_EDGES["c5"])
    candidate_modulated, _ = simulate_sequence(
        candidate,
        GRAPH_EDGES["c5"],
        with_modulation=True,
    )
    report = {
        "fallback": {
            "ideal_p_mis": fallback_ideal,
            "modulated_p_mis": fallback_modulated,
        },
        "robust_candidate": {
            "dense_nominal_p_mis": dense_scores[0],
            "dense_worst_case_p_mis": min(dense_scores),
            "dense_ensemble_p_mis": dense_scores,
            "ideal_p_mis": candidate_ideal,
            "modulated_p_mis": candidate_modulated,
            "improves_modulated_score": (
                candidate_modulated > fallback_modulated
            ),
        },
        "parameters": parameters,
    }
    save_sequence(candidate, output / "sequence_c5_robust.json")
    save_json(parameters, output / "parameters_c5_robust.json")
    save_json(report, output / "scores_c5_robust.json")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Optimize a modulation-safe, ensemble-robust C5 pulse."
    )
    parser.add_argument("--output", default="results/challenge02")
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    report = run(args.output, args.quick)
    fallback = report["fallback"]
    robust = report["robust_candidate"]
    print(
        f"fallback modulated={fallback['modulated_p_mis']:.8f}",
        f"robust modulated={robust['modulated_p_mis']:.8f}",
        f"robust worst_dense={robust['dense_worst_case_p_mis']:.8f}",
    )


if __name__ == "__main__":
    main()
