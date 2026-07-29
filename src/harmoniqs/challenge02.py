from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from pulser import AnalogDevice, Pulse, Register, Sequence
from pulser.waveforms import (
    CompositeWaveform,
    ConstantWaveform,
    InterpolatedWaveform,
    RampWaveform,
)
from pulser_simulation import QutipEmulator
from scipy.linalg import eigh
from scipy.optimize import differential_evolution, minimize

from .common import (
    C6_OVER_HBAR,
    deterministic_counts,
    mis_bitstrings,
    mis_probability,
    probabilities_from_state,
    register_with_layout,
    save_json,
    save_sequence,
    validate_sequence,
)

RISE_NS = 252
REFERENCE = {
    "omega": 2 * np.pi,
    "detuning_initial": -4 * np.pi,
    "detuning_final": 4 * np.pi,
    "power": 1.0,
    "duration_ns": 4000,
    "spacing_um": 5.5,
}
GRAPH_EDGES = {
    "star": ((0, 1), (0, 2), (0, 3)),
    "c5": ((0, 1), (1, 2), (2, 3), (3, 4), (0, 4)),
}


def register_coordinates(graph: str, spacing_um: float) -> list[tuple[float, float]]:
    if graph == "star":
        leaves = [
            (
                spacing_um * np.cos(2 * np.pi * k / 3),
                spacing_um * np.sin(2 * np.pi * k / 3),
            )
            for k in range(3)
        ]
        return [(0.0, 0.0), *leaves]
    if graph == "c5":
        radius = spacing_um / (2 * np.sin(np.pi / 5))
        return [
            (
                radius * np.cos(np.pi / 2 + 2 * np.pi * k / 5),
                radius * np.sin(np.pi / 2 + 2 * np.pi * k / 5),
            )
            for k in range(5)
        ]
    raise ValueError(f"Unknown graph: {graph}")


def make_register(graph: str, spacing_um: float) -> Register:
    coordinates = register_coordinates(graph, spacing_um)
    return register_with_layout(
        coordinates,
        [f"q{k}" for k in range(len(coordinates))],
    )


def make_sequence(graph: str, parameters: dict[str, float | int]) -> Sequence:
    duration_ns = int(parameters["duration_ns"])
    hold_ns = duration_ns - 2 * RISE_NS
    omega = float(parameters["omega"])
    detuning_initial = float(parameters["detuning_initial"])
    detuning_final = float(parameters["detuning_final"])
    power = float(parameters["power"])
    spacing = float(parameters["spacing_um"])

    amplitude = CompositeWaveform(
        RampWaveform(RISE_NS, 0.0, omega),
        ConstantWaveform(hold_ns, omega),
        RampWaveform(RISE_NS, omega, 0.0),
    )
    times = np.linspace(0.0, 1.0, 17)
    detuning_values = detuning_initial + (
        detuning_final - detuning_initial
    ) * times**power
    detuning = CompositeWaveform(
        ConstantWaveform(RISE_NS, detuning_initial),
        InterpolatedWaveform(hold_ns, detuning_values, times=times),
        ConstantWaveform(RISE_NS, detuning_final),
    )
    sequence = Sequence(make_register(graph, spacing), AnalogDevice)
    sequence.declare_channel("rydberg", "rydberg_global")
    sequence.add(Pulse(amplitude, detuning, 0.0), "rydberg")
    validate_sequence(sequence)
    return sequence


def _operators(
    coordinates: list[tuple[float, float]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(coordinates)
    dimension = 2**n
    drive = np.zeros((dimension, dimension), dtype=float)
    excitations = np.zeros(dimension, dtype=float)
    interactions = np.zeros(dimension, dtype=float)
    for basis in range(dimension):
        bits = format(basis, f"0{n}b")
        excitations[basis] = bits.count("1")
        for qubit in range(n):
            drive[basis, basis ^ (1 << (n - 1 - qubit))] += 0.5
        for i in range(n):
            if bits[i] != "1":
                continue
            for j in range(i + 1, n):
                if bits[j] != "1":
                    continue
                distance = np.linalg.norm(
                    np.asarray(coordinates[i]) - np.asarray(coordinates[j])
                )
                interactions[basis] += C6_OVER_HBAR / distance**6
    return drive, excitations, interactions


def control_at(time_us: float, parameters: dict[str, float | int]) -> tuple[float, float]:
    duration_us = float(parameters["duration_ns"]) / 1000
    rise_us = RISE_NS / 1000
    hold_us = duration_us - 2 * rise_us
    omega_max = float(parameters["omega"])
    initial = float(parameters["detuning_initial"])
    final = float(parameters["detuning_final"])
    if time_us < rise_us:
        return omega_max * time_us / rise_us, initial
    if time_us <= rise_us + hold_us:
        x = (time_us - rise_us) / hold_us
        return omega_max, initial + (final - initial) * x ** float(parameters["power"])
    return omega_max * (duration_us - time_us) / rise_us, final


def evolve_dense(
    graph: str, parameters: dict[str, float | int], steps: int = 72
) -> np.ndarray:
    coordinates = register_coordinates(graph, float(parameters["spacing_um"]))
    drive, excitations, interactions = _operators(coordinates)
    duration_us = float(parameters["duration_ns"]) / 1000
    dt = duration_us / steps
    state = np.zeros(drive.shape[0], dtype=complex)
    state[0] = 1.0
    for step in range(steps):
        omega, detuning = control_at((step + 0.5) * dt, parameters)
        hamiltonian = omega * drive + np.diag(interactions - detuning * excitations)
        eigenvalues, eigenvectors = eigh(hamiltonian)
        state = eigenvectors @ (
            np.exp(-1j * eigenvalues * dt) * (eigenvectors.conj().T @ state)
        )
    return state


def simulate_sequence(
    sequence: Sequence,
    edges: tuple[tuple[int, int], ...],
    *,
    with_modulation: bool = False,
) -> tuple[float, dict[str, float]]:
    final = (
        QutipEmulator.from_sequence(
            sequence,
            evaluation_times=0.01,
            with_modulation=with_modulation,
        )
        .run()
        .get_final_state()
        .full()
        .reshape(-1)
    )
    probabilities = probabilities_from_state(final, pulser_basis=True)
    return mis_probability(probabilities, edges), probabilities


def optimize_parameters(
    graph: str, *, quick: bool = False, seed: int = 73
) -> tuple[dict[str, float | int], float]:
    edges = GRAPH_EDGES[graph]
    bounds = [
        (4.5, 10.0),
        (-24.0, -5.0),
        (6.0, 24.0),
        (0.35, 2.5),
        (5.02, 6.8),
    ]

    def unpack(values: np.ndarray) -> dict[str, float | int]:
        return {
            "omega": float(values[0]),
            "detuning_initial": float(values[1]),
            "detuning_final": float(values[2]),
            "power": float(values[3]),
            "duration_ns": 6000,
            "spacing_um": float(values[4]),
        }

    targets = mis_bitstrings(
        4 if graph == "star" else 5,
        edges,
    )

    def objective(values: np.ndarray) -> float:
        state = evolve_dense(graph, unpack(values), steps=44 if quick else 72)
        probabilities = probabilities_from_state(state)
        return -sum(probabilities[bits] for bits in targets)

    initial = np.array(
        [
            REFERENCE["omega"],
            REFERENCE["detuning_initial"],
            REFERENCE["detuning_final"],
            1.0,
            5.5,
        ]
    )
    result = differential_evolution(
        objective,
        bounds,
        x0=initial,
        seed=seed,
        maxiter=8 if quick else 24,
        popsize=5 if quick else 7,
        polish=False,
        workers=1,
    )
    polished = minimize(
        objective,
        result.x,
        method="Powell",
        bounds=bounds,
        options={
            "maxiter": 45 if quick else 120,
            "xtol": 2e-4,
            "ftol": 1e-9,
        },
    )
    clipped = np.array(
        [np.clip(value, low, high) for value, (low, high) in zip(polished.x, bounds)]
    )
    parameters = unpack(clipped)
    score = -objective(clipped)
    return parameters, float(score)


def run(
    output_dir: str | Path = "results/challenge02",
    graphs: tuple[str, ...] = ("star", "c5"),
    quick: bool = False,
) -> dict:
    output = Path(output_dir)
    report: dict[str, object] = {"challenge": 2, "graphs": {}}
    for graph_index, graph in enumerate(graphs):
        edges = GRAPH_EDGES[graph]
        baseline_sequence = make_sequence(graph, REFERENCE)
        baseline_score, baseline_probabilities = simulate_sequence(
            baseline_sequence, edges
        )
        parameters, dense_score = optimize_parameters(
            graph, quick=quick, seed=73 + graph_index
        )
        optimized_sequence = make_sequence(graph, parameters)
        optimized_score, optimized_probabilities = simulate_sequence(
            optimized_sequence, edges
        )
        modulated_score, modulated_probabilities = simulate_sequence(
            optimized_sequence, edges, with_modulation=True
        )
        save_sequence(optimized_sequence, output / f"sequence_{graph}.json")
        save_json(parameters, output / f"parameters_{graph}.json")
        report["graphs"][graph] = {
            "maximum_independent_sets": sorted(
                mis_bitstrings(len(optimized_sequence.register.qubit_ids), edges)
            ),
            "baseline_p_mis": baseline_score,
            "optimized_dense_p_mis": dense_score,
            "optimized_p_mis": optimized_score,
            "optimized_modulated_p_mis": modulated_score,
            "strictly_beats_baseline": optimized_score > baseline_score,
            "baseline_counts_500": dict(
                deterministic_counts(baseline_probabilities)
            ),
            "optimized_counts_500": dict(
                deterministic_counts(optimized_probabilities)
            ),
            "optimized_modulated_counts_500": dict(
                deterministic_counts(modulated_probabilities)
            ),
            "parameters": parameters,
        }
    save_json(report, output / "scores.json")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Optimize Challenge 02.")
    parser.add_argument("--output", default="results/challenge02")
    parser.add_argument(
        "--graph", choices=("star", "c5", "both"), default="both"
    )
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    graphs = ("star", "c5") if args.graph == "both" else (args.graph,)
    report = run(args.output, graphs, args.quick)
    for graph, values in report["graphs"].items():
        print(
            graph,
            f"baseline={values['baseline_p_mis']:.6f}",
            f"optimized={values['optimized_p_mis']:.6f}",
            f"modulated={values['optimized_modulated_p_mis']:.6f}",
        )


if __name__ == "__main__":
    main()
