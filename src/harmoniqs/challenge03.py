from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence as TypingSequence

import numpy as np
from pulser import AnalogDevice, Pulse, Sequence
from pulser.sampler import sample
from pulser.waveforms import (
    CompositeWaveform,
    ConstantWaveform,
    InterpolatedWaveform,
    RampWaveform,
)
from scipy.interpolate import PchipInterpolator
from scipy.optimize import differential_evolution, minimize
from scipy.sparse import csr_matrix

from .common import (
    C6_OVER_HBAR,
    deterministic_counts,
    register_with_layout,
    save_json,
    save_sequence,
    validate_sequence,
)

PAPER_REPOSITORY = "https://github.com/alejomonbar/Benchmarking-neutral-atom-QPUs"
PAPER_TARGETS = {11: 0.90661639, 13: 0.90817732, 17: 0.87004419}
PAPER_VALID_FRACTIONS = {11: 0.65266667, 13: 0.66133333, 17: 0.515}
PAPER_500_SHOT_REPLICATES = {
    11: (
        {"case": 2, "valid_count": 276, "valid_size_sum": 1011},
        {"case": 3, "valid_count": 356, "valid_size_sum": 1278},
    ),
    13: (
        {"case": 2, "valid_count": 294, "valid_size_sum": 1088},
        {"case": 3, "valid_count": 363, "valid_size_sum": 1294},
    ),
    17: (
        {"case": 2, "valid_count": 240, "valid_size_sum": 1249},
        {"case": 3, "valid_count": 275, "valid_size_sum": 1440},
    ),
}
GRID_SPACING_UM = 5.0
KNOT_TIMES = np.linspace(0.0, 1.0, 9)

# Exact Data/Problems/{11,13}.json instances from the paper repository. Vertex
# order is significant because it defines the measurement bitstring order.
_POSITIONS_11 = (
    (3, 1), (2, 0), (2, 2), (0, 2), (0, 0), (3, 3),
    (2, 3), (0, 1), (0, 3), (1, 0), (1, 1),
)
_EDGES_11 = (
    (0, 1), (0, 2), (1, 9), (1, 10), (2, 5), (2, 6), (2, 10),
    (3, 7), (3, 8), (3, 10), (4, 7), (4, 9), (4, 10), (5, 6),
    (7, 9), (7, 10), (9, 10),
)
_POSITIONS_13 = _POSITIONS_11 + ((2, 1), (3, 2))
_EDGES_13 = (
    (0, 1), (0, 2), (0, 11), (0, 12), (1, 9), (1, 10), (1, 11),
    (2, 5), (2, 6), (2, 10), (2, 11), (2, 12), (3, 7), (3, 8),
    (3, 10), (4, 7), (4, 9), (4, 10), (5, 6), (5, 12), (6, 12),
    (7, 9), (7, 10), (9, 10), (9, 11), (10, 11), (11, 12),
)
_POSITIONS_17 = (
    (0, 4), (2, 0), (4, 0), (4, 4), (0, 2), (3, 4), (2, 1),
    (1, 1), (1, 2), (1, 0), (0, 0), (3, 1), (3, 3), (1, 4),
    (0, 3), (4, 3), (1, 3),
)
_EDGES_17 = (
    (0, 13), (0, 14), (0, 16), (1, 6), (1, 7), (1, 9), (1, 11),
    (2, 11), (3, 5), (3, 12), (3, 15), (4, 7), (4, 8), (4, 14),
    (4, 16), (5, 12), (5, 15), (6, 7), (6, 8), (6, 9), (6, 11),
    (7, 8), (7, 9), (7, 10), (8, 14), (8, 16), (9, 10), (12, 15),
    (13, 14), (13, 16), (14, 16),
)


@dataclass(frozen=True)
class BenchmarkInstance:
    n_vertices: int
    grid_side: int
    grid_positions: tuple[tuple[int, int], ...]
    edges: tuple[tuple[int, int], ...]
    paper_target: float
    paper_valid_fraction: float

    def coordinates(self, spacing_um: float = GRID_SPACING_UM) -> list[tuple[float, float]]:
        shift = (self.grid_side - 1) * spacing_um / 2
        return [
            (spacing_um * x - shift, spacing_um * y - shift)
            for x, y in self.grid_positions
        ]


INSTANCES = {
    11: BenchmarkInstance(
        11, 4, _POSITIONS_11, _EDGES_11, PAPER_TARGETS[11], PAPER_VALID_FRACTIONS[11]
    ),
    13: BenchmarkInstance(
        13, 4, _POSITIONS_13, _EDGES_13, PAPER_TARGETS[13], PAPER_VALID_FRACTIONS[13]
    ),
    17: BenchmarkInstance(
        17, 5, _POSITIONS_17, _EDGES_17, PAPER_TARGETS[17], PAPER_VALID_FRACTIONS[17]
    ),
}


def load_instance(n_vertices: int) -> BenchmarkInstance:
    try:
        return INSTANCES[n_vertices]
    except KeyError as exc:
        raise ValueError(
            "Challenge 3 supports the exact N=11, N=13, and N=17 instances."
        ) from exc


def is_independent(bits: str, edges: TypingSequence[tuple[int, int]]) -> bool:
    return all(not (bits[i] == bits[j] == "1") for i, j in edges)


def qubo_cost(bits: str, edges: TypingSequence[tuple[int, int]]) -> int:
    violations = sum(bits[i] == bits[j] == "1" for i, j in edges)
    return -bits.count("1") + 2 * violations


def certify_instance(instance: BenchmarkInstance) -> dict[str, object]:
    valid = [
        format(value, f"0{instance.n_vertices}b")
        for value in range(2**instance.n_vertices)
        if is_independent(format(value, f"0{instance.n_vertices}b"), instance.edges)
    ]
    alpha = max(bits.count("1") for bits in valid)
    optima = [bits for bits in valid if bits.count("1") == alpha]
    return {
        "n_vertices": instance.n_vertices,
        "alpha": alpha,
        "minimum_qubo_cost": -alpha,
        "number_of_maximum_independent_sets": len(optima),
        "maximum_independent_sets": optima,
        "edge_count": len(instance.edges),
    }


def score_distribution(
    probabilities: Mapping[str, float],
    instance: BenchmarkInstance,
) -> dict[str, float]:
    certificate = certify_instance(instance)
    alpha = int(certificate["alpha"])
    valid_probability = 0.0
    valid_size = 0.0
    optimal_probability = 0.0
    expected_cost = 0.0
    for bits, probability in probabilities.items():
        expected_cost += probability * qubo_cost(bits, instance.edges)
        if is_independent(bits, instance.edges):
            valid_probability += probability
            valid_size += probability * bits.count("1")
            if bits.count("1") == alpha:
                optimal_probability += probability
    ratio = valid_size / (alpha * valid_probability) if valid_probability else 0.0
    return {
        "valid_probability": float(valid_probability),
        "valid_approximation_ratio": float(ratio),
        "optimal_probability": float(optimal_probability),
        "expected_qubo_cost": float(expected_cost),
    }


def score_counts(
    counts: Mapping[str, int],
    instance: BenchmarkInstance,
) -> dict[str, float | int]:
    shots = int(sum(counts.values()))
    probabilities = {bits: count / shots for bits, count in counts.items()}
    return {"shots": shots, **score_distribution(probabilities, instance)}


@dataclass(frozen=True)
class SparseRydbergModel:
    coordinates: np.ndarray
    drive: csr_matrix
    excitations: np.ndarray
    interactions: np.ndarray
    drive_pairs: tuple[tuple[np.ndarray, np.ndarray], ...]

    @classmethod
    def from_coordinates(
        cls, coordinates: TypingSequence[tuple[float, float]]
    ) -> "SparseRydbergModel":
        coordinates_array = np.asarray(coordinates, dtype=float)
        n = len(coordinates_array)
        dimension = 2**n
        basis = np.arange(dimension, dtype=np.int64)
        rows = np.repeat(basis, n)
        cols = np.concatenate(
            [basis ^ (1 << (n - 1 - qubit)) for qubit in range(n)]
        ).reshape(n, dimension).T.reshape(-1)
        drive = csr_matrix(
            (np.full(rows.size, 0.5), (rows, cols)),
            shape=(dimension, dimension),
        )
        bit_masks = 1 << np.arange(n - 1, -1, -1)
        occupied = (basis[:, None] & bit_masks[None, :]) != 0
        excitations = occupied.sum(axis=1).astype(float)
        interactions = np.zeros(dimension, dtype=float)
        for i in range(n):
            for j in range(i + 1, n):
                distance = np.linalg.norm(coordinates_array[i] - coordinates_array[j])
                interactions += (
                    C6_OVER_HBAR / distance**6 * (occupied[:, i] & occupied[:, j])
                )
        drive_pairs = tuple(
            (
                basis[(basis & bit_mask) == 0],
                basis[(basis & bit_mask) != 0],
            )
            for bit_mask in bit_masks
        )
        return cls(coordinates_array, drive, excitations, interactions, drive_pairs)


def controls_at(
    times: np.ndarray, parameters: Mapping[str, object]
) -> tuple[np.ndarray, np.ndarray]:
    duration_us = float(parameters["duration_ns"]) / 1000.0
    normalized = np.clip(times / duration_us, 0.0, 1.0)
    knot_times = np.asarray(parameters.get("knot_times", KNOT_TIMES), dtype=float)
    amplitude_knots = np.asarray(parameters["amplitude_knots"], dtype=float)
    detuning_knots = np.asarray(parameters["detuning_knots"], dtype=float)
    if parameters.get("interpolator") == "interp1d":
        amplitude = np.interp(normalized, knot_times, amplitude_knots)
        detuning = np.interp(normalized, knot_times, detuning_knots)
    else:
        amplitude = PchipInterpolator(knot_times, amplitude_knots)(normalized)
        detuning = PchipInterpolator(knot_times, detuning_knots)(normalized)
    return np.maximum(amplitude, 0.0), detuning


def evolve_sparse(
    model: SparseRydbergModel,
    parameters: Mapping[str, object],
    *,
    steps: int = 120,
    omega_scale: float = 1.0,
    detuning_offset: float = 0.0,
) -> np.ndarray:
    duration_us = float(parameters["duration_ns"]) / 1000.0
    dt = duration_us / steps
    midpoints = (np.arange(steps) + 0.5) * dt
    amplitudes, detunings = controls_at(midpoints, parameters)
    state = np.zeros(model.drive.shape[0], dtype=complex)
    state[0] = 1.0
    for amplitude, detuning in zip(amplitudes, detunings):
        state = _split_step(
            model,
            state,
            omega_scale * amplitude,
            detuning + detuning_offset,
            dt,
        )
    return state


def _split_step(
    model: SparseRydbergModel,
    state: np.ndarray,
    amplitude: float,
    detuning: float,
    dt_us: float,
) -> np.ndarray:
    """Apply one second-order diagonal/global-drive split step."""
    diagonal = model.interactions - detuning * model.excitations
    half_phase = np.exp(-0.5j * dt_us * diagonal)
    state *= half_phase
    angle = amplitude * dt_us / 2.0
    cosine, sine = np.cos(angle), -1j * np.sin(angle)
    for zero, one in model.drive_pairs:
        ground = state[zero].copy()
        excited = state[one].copy()
        state[zero] = cosine * ground + sine * excited
        state[one] = sine * ground + cosine * excited
    state *= half_phase
    return state


def modulated_controls(
    sequence: Sequence,
    *,
    bin_ns: int = 16,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sample Pulser's channel model and average it into integration bins."""
    nested = sample(sequence, modulation=True).to_nested_dict()
    channel = nested["Global"]["ground-rydberg"]
    amplitude = np.asarray(channel["amp"], dtype=float)
    detuning = np.asarray(channel["det"], dtype=float)
    amplitudes: list[float] = []
    detunings: list[float] = []
    durations: list[float] = []
    for start in range(0, len(amplitude), bin_ns):
        stop = min(start + bin_ns, len(amplitude))
        amplitudes.append(float(np.mean(amplitude[start:stop])))
        detunings.append(float(np.mean(detuning[start:stop])))
        durations.append((stop - start) / 1000.0)
    return np.asarray(amplitudes), np.asarray(detunings), np.asarray(durations)


def evolve_modulated(
    model: SparseRydbergModel,
    sequence: Sequence,
    *,
    bin_ns: int = 16,
) -> np.ndarray:
    amplitudes, detunings, durations = modulated_controls(sequence, bin_ns=bin_ns)
    state = np.zeros(model.drive.shape[0], dtype=complex)
    state[0] = 1.0
    for amplitude, detuning, duration in zip(amplitudes, detunings, durations):
        state = _split_step(model, state, amplitude, detuning, duration)
    return state


def probabilities_from_state(state: np.ndarray, n_vertices: int) -> dict[str, float]:
    return {
        format(index, f"0{n_vertices}b"): float(abs(amplitude) ** 2)
        for index, amplitude in enumerate(state)
    }


def make_register(instance: BenchmarkInstance, spacing_um: float = GRID_SPACING_UM):
    return register_with_layout(
        instance.coordinates(spacing_um),
        [f"q{index}" for index in range(instance.n_vertices)],
    )


def make_sequence(
    instance: BenchmarkInstance,
    parameters: Mapping[str, object],
) -> Sequence:
    duration_ns = int(parameters["duration_ns"])
    if parameters.get("interpolator") == "interp1d":
        knot_times = np.asarray(parameters["knot_times"], dtype=float)
        boundaries = np.rint(knot_times * duration_ns).astype(int)
        durations = np.diff(boundaries)
        amplitude_values = np.asarray(parameters["amplitude_knots"], dtype=float)
        detuning_values = np.asarray(parameters["detuning_knots"], dtype=float)
        amplitude = CompositeWaveform(
            *[
                RampWaveform(int(segment), amplitude_values[i], amplitude_values[i + 1])
                for i, segment in enumerate(durations)
            ]
        )
        detuning = CompositeWaveform(
            *[
                RampWaveform(int(segment), detuning_values[i], detuning_values[i + 1])
                for i, segment in enumerate(durations)
            ]
        )
    else:
        knot_times = np.asarray(parameters.get("knot_times", KNOT_TIMES), dtype=float)
        amplitude = InterpolatedWaveform(
            duration_ns,
            np.asarray(parameters["amplitude_knots"], dtype=float),
            times=knot_times,
        )
        detuning = InterpolatedWaveform(
            duration_ns,
            np.asarray(parameters["detuning_knots"], dtype=float),
            times=knot_times,
        )
    sequence = Sequence(
        make_register(instance, float(parameters.get("spacing_um", GRID_SPACING_UM))),
        AnalogDevice,
    )
    sequence.declare_channel("rydberg", "rydberg_global")
    sequence.add(Pulse(amplitude, detuning, 0.0), "rydberg")
    validate_sequence(sequence)
    return sequence


def paper_reference_parameters() -> dict[str, object]:
    # Exact 4 us paper schedule: 1420 ns rise, 1160 ns sweep, 1420 ns fall.
    return {
        "duration_ns": 4000,
        "spacing_um": GRID_SPACING_UM,
        "knot_times": [0.0, 1420 / 4000, 2580 / 4000, 1.0],
        "interpolator": "interp1d",
        "amplitude_knots": [0.0, 5.5, 5.5, 0.0],
        "detuning_knots": [-5.5, -5.5, 5.5, 5.5],
    }


def scaled_reference_parameters() -> dict[str, object]:
    """Paper QAA shape stretched from 4 µs to the 6 µs device limit."""
    return {
        "duration_ns": 6000,
        "spacing_um": GRID_SPACING_UM,
        "knot_times": [0.0, 0.355, 0.645, 1.0],
        "interpolator": "interp1d",
        "amplitude_knots": [0.0, 5.5, 5.5, 0.0],
        "detuning_knots": [-5.5, -5.5, 5.5, 5.5],
    }


def smooth_parameters(
    values: TypingSequence[float],
    *,
    duration_ns: int = 6000,
) -> dict[str, object]:
    """Build a smooth trapezoid with a shaped monotone detuning sweep."""
    omega, initial, final, rise_fraction, sweep_power = map(float, values)
    times = KNOT_TIMES
    amplitude = np.ones_like(times)
    rising = times < rise_fraction
    falling = times > 1.0 - rise_fraction
    rise_x = np.clip(times[rising] / rise_fraction, 0.0, 1.0)
    fall_x = np.clip((1.0 - times[falling]) / rise_fraction, 0.0, 1.0)
    amplitude[rising] = rise_x * rise_x * (3.0 - 2.0 * rise_x)
    amplitude[falling] = fall_x * fall_x * (3.0 - 2.0 * fall_x)
    sweep_x = np.clip(
        (times - rise_fraction) / (1.0 - 2.0 * rise_fraction), 0.0, 1.0
    )
    left = np.power(sweep_x, sweep_power)
    right = np.power(1.0 - sweep_x, sweep_power)
    shaped = np.divide(left, left + right, out=np.zeros_like(left), where=left + right > 0)
    return {
        "duration_ns": duration_ns,
        "spacing_um": GRID_SPACING_UM,
        "amplitude_knots": (omega * amplitude).tolist(),
        "detuning_knots": (initial + (final - initial) * shaped).tolist(),
    }


def unpack_parameters(values: np.ndarray) -> dict[str, object]:
    omega, initial, final = values[:3]
    amplitude_logits = values[3:6]
    amplitude_interior = omega / (1.0 + np.exp(-amplitude_logits))
    amplitude_knots = np.concatenate(
        ([0.0], amplitude_interior, [omega], amplitude_interior[::-1], [0.0])
    )
    detuning_logits = np.append(values[6:], 0.0)
    weights = np.exp(detuning_logits - np.max(detuning_logits))
    weights /= weights.sum()
    detuning_knots = initial + (final - initial) * np.concatenate(
        ([0.0], np.cumsum(weights))
    )
    return {
        "duration_ns": 6000,
        "spacing_um": GRID_SPACING_UM,
        "amplitude_knots": amplitude_knots.tolist(),
        "detuning_knots": detuning_knots.tolist(),
    }


def optimize(
    n_vertices: int = 11,
    *,
    quick: bool = False,
    seed: int = 311,
) -> tuple[dict[str, object], dict[str, float]]:
    instances = [load_instance(n) for n in (11, 13)]
    models = [
        SparseRydbergModel.from_coordinates(instance.coordinates())
        for instance in instances
    ]
    bounds = [
        (4.0, 8.0),
        (-10.0, -3.0),
        (3.0, 12.0),
        (0.18, 0.40),
        (0.5, 2.5),
    ]
    initial = np.array([5.5, -5.5, 5.5, 0.30, 1.0])
    search_steps = 64 if quick else 128

    def objective(values: np.ndarray) -> float:
        parameters = smooth_parameters(values)
        scores = []
        for instance, model in zip(instances, models):
            state = evolve_sparse(model, parameters, steps=search_steps)
            scores.append(
                score_distribution(
                    probabilities_from_state(state, instance.n_vertices), instance
                )
            )
        ratios = np.array([score["valid_approximation_ratio"] for score in scores])
        valid = np.array([score["valid_probability"] for score in scores])
        optimum = np.array([score["optimal_probability"] for score in scores])
        quality = 0.45 * ratios.mean() + 0.20 * ratios.min()
        quality += 0.20 * valid.mean() + 0.15 * optimum.mean()
        return -float(quality)

    result = differential_evolution(
        objective,
        bounds,
        x0=initial,
        seed=seed,
        maxiter=3 if quick else 12,
        popsize=3 if quick else 5,
        workers=1,
        polish=False,
        updating="immediate",
    )
    polished = minimize(
        objective,
        result.x,
        method="Powell",
        bounds=bounds,
        options={"maxiter": 20 if quick else 70, "ftol": 1e-8},
    )
    parameters = smooth_parameters(polished.x)
    selected_instance = load_instance(n_vertices)
    selected_model = SparseRydbergModel.from_coordinates(selected_instance.coordinates())
    final_state = evolve_sparse(selected_model, parameters, steps=384)
    score = score_distribution(
        probabilities_from_state(final_state, n_vertices), selected_instance
    )
    return parameters, score


def optimize_hardware_retry(
    *,
    seed: int = 317,
    quick: bool = False,
) -> tuple[dict[str, object], dict[str, float]]:
    """Optimize a shorter N=17 pulse after the first hardware calibration run."""
    instance = load_instance(17)
    model = SparseRydbergModel.from_coordinates(instance.coordinates())
    bounds = [
        (3.5, 6.5),
        (-8.0, -3.0),
        (4.0, 10.0),
        (0.22, 0.40),
        (0.45, 1.8),
    ]
    initial = np.array([4.0000045, -4.648809, 6.590995, 0.37, 0.5])
    steps = 72 if quick else 96

    def objective(values: np.ndarray) -> float:
        parameters = smooth_parameters(values, duration_ns=4000)
        state = evolve_sparse(model, parameters, steps=steps)
        score = score_distribution(
            probabilities_from_state(state, instance.n_vertices), instance
        )
        validity_shortfall = max(0.0, 0.65 - score["valid_probability"])
        quality = score["valid_approximation_ratio"]
        quality += 0.02 * score["optimal_probability"]
        quality -= 2.0 * validity_shortfall**2
        return -float(quality)

    result = differential_evolution(
        objective,
        bounds,
        x0=initial,
        seed=seed,
        maxiter=1 if quick else 3,
        popsize=2,
        workers=1,
        polish=False,
        updating="immediate",
    )
    values = result.x if objective(result.x) < objective(initial) else initial
    parameters = smooth_parameters(values, duration_ns=4000)
    modulated, _ = evaluate_modulated(instance, parameters, bin_ns=24)
    return parameters, modulated


def run_hardware_retry(
    output_dir: str | Path = "results/challenge03",
    *,
    quick: bool = False,
) -> dict[str, object]:
    output = Path(output_dir)
    instance = load_instance(17)
    parameters, search_score = optimize_hardware_retry(quick=quick)
    sequence = make_sequence(instance, parameters)
    ideal, ideal_500 = evaluate(instance, parameters, steps=192)
    modulated, modulated_500 = evaluate_modulated(instance, parameters, bin_ns=24)
    model = SparseRydbergModel.from_coordinates(instance.coordinates())
    robust = []
    for label, omega_scale, detuning_offset in (
        ("omega_-2pct", 0.98, 0.0),
        ("omega_+2pct", 1.02, 0.0),
        ("detuning_-0.3", 1.0, -0.3),
        ("detuning_+0.3", 1.0, 0.3),
    ):
        state = evolve_sparse(
            model,
            parameters,
            steps=128,
            omega_scale=omega_scale,
            detuning_offset=detuning_offset,
        )
        robust.append(
            {
                "variation": label,
                "score": score_distribution(
                    probabilities_from_state(state, instance.n_vertices), instance
                ),
            }
        )
    report = {
        "purpose": "N=17 hardware retry after calibration run",
        "hardware_submitted": False,
        "duration_ns": 4000,
        "paper_target": instance.paper_target,
        "first_hardware_score": 0.8452380952380952,
        "search_modulated_score": search_score,
        "ideal_exact": ideal,
        "ideal_counts_500": ideal_500,
        "modulated_exact": modulated,
        "modulated_counts_500": modulated_500,
        "robust_ensemble": robust,
        "robust_worst_valid_approximation_ratio": min(
            entry["score"]["valid_approximation_ratio"] for entry in robust
        ),
        "ready_for_hardware": (
            modulated["valid_approximation_ratio"] > 0.96
            and modulated["valid_probability"] > 0.60
            and min(
                entry["score"]["valid_approximation_ratio"] for entry in robust
            )
            > 0.94
        ),
    }
    save_sequence(sequence, output / "sequence_n17_retry.json")
    save_json(parameters, output / "parameters_n17_retry.json")
    save_json(report, output / "scores_n17_retry.json")
    return report


def evaluate(
    instance: BenchmarkInstance,
    parameters: Mapping[str, object],
    *,
    steps: int = 192,
) -> tuple[dict[str, float], dict[str, float]]:
    model = SparseRydbergModel.from_coordinates(instance.coordinates())
    state = evolve_sparse(model, parameters, steps=steps)
    probabilities = probabilities_from_state(state, instance.n_vertices)
    exact = score_distribution(probabilities, instance)
    counts = deterministic_counts(probabilities)
    shot_score = score_counts(counts, instance)
    return exact, shot_score


def evaluate_modulated(
    instance: BenchmarkInstance,
    parameters: Mapping[str, object],
    *,
    bin_ns: int = 16,
) -> tuple[dict[str, float], dict[str, float]]:
    model = SparseRydbergModel.from_coordinates(instance.coordinates())
    sequence = make_sequence(instance, parameters)
    state = evolve_modulated(model, sequence, bin_ns=bin_ns)
    probabilities = probabilities_from_state(state, instance.n_vertices)
    return (
        score_distribution(probabilities, instance),
        score_counts(deterministic_counts(probabilities), instance),
    )


def run(
    output_dir: str | Path = "results/challenge03",
    *,
    quick: bool = False,
) -> dict[str, object]:
    output = Path(output_dir)
    parameters, optimized_n11 = optimize(11, quick=quick)
    summary: dict[str, object] = {
        "challenge": 3,
        "source": PAPER_REPOSITORY,
        "hardware_submitted": False,
        "instances": {},
    }
    save_json(
        {
            "paper": "arXiv:2511.22967",
            "repository": PAPER_REPOSITORY,
            "instance_paths": [
                "Data/Problems/11.json",
                "Data/Problems/13.json",
                "Data/Problems/17.json",
            ],
            "generation": "4x4/5x5 diagonal-connected grids, NumPy seed 123, random dropout",
            "paper_scoring": "mean selected-set size over valid samples divided by alpha(G)",
            "paper_replicates_used": [2, 3],
            "hardware_jobs_submitted": False,
        },
        output / "provenance.json",
    )
    for n_vertices in (11, 13, 17):
        instance = load_instance(n_vertices)
        certificate = certify_instance(instance)
        alpha = int(certificate["alpha"])
        final_steps = 192 if n_vertices == 17 else 384
        robust_steps = 96 if n_vertices == 17 else 192
        paper_replicates = [
            {
                **replicate,
                "shots": 500,
                "valid_approximation_ratio": (
                    replicate["valid_size_sum"] / (alpha * replicate["valid_count"])
                ),
            }
            for replicate in PAPER_500_SHOT_REPLICATES[n_vertices]
        ]
        baseline_exact, baseline_500 = evaluate(
            instance, paper_reference_parameters(), steps=final_steps
        )
        optimized_exact, optimized_500 = evaluate(
            instance, parameters, steps=final_steps
        )
        modulated_exact, modulated_500 = evaluate_modulated(
            instance,
            parameters,
            bin_ns=32 if n_vertices == 17 else 16,
        )
        model = SparseRydbergModel.from_coordinates(instance.coordinates())
        convergence = {}
        convergence_steps = (48, 96, 192) if n_vertices == 17 else (96, 192, 384)
        for steps in convergence_steps:
            state = evolve_sparse(model, parameters, steps=steps)
            convergence[str(steps)] = score_distribution(
                probabilities_from_state(state, n_vertices), instance
            )
        robust_scores: list[dict[str, object]] = []
        for label, omega_scale, detuning_offset in (
            ("omega_-2pct", 0.98, 0.0),
            ("omega_+2pct", 1.02, 0.0),
            ("detuning_-0.3", 1.0, -0.3),
            ("detuning_+0.3", 1.0, 0.3),
        ):
            state = evolve_sparse(
                model,
                parameters,
                steps=robust_steps,
                omega_scale=omega_scale,
                detuning_offset=detuning_offset,
            )
            robust_scores.append(
                {
                    "variation": label,
                    "score": score_distribution(
                        probabilities_from_state(state, n_vertices), instance
                    ),
                }
            )
        for spacing_offset in (-0.02, 0.02):
            perturbed_model = SparseRydbergModel.from_coordinates(
                instance.coordinates(GRID_SPACING_UM + spacing_offset)
            )
            state = evolve_sparse(perturbed_model, parameters, steps=robust_steps)
            robust_scores.append(
                {
                    "variation": f"spacing_{spacing_offset:+.2f}um",
                    "score": score_distribution(
                        probabilities_from_state(state, n_vertices), instance
                    ),
                }
            )
        sequence = make_sequence(instance, parameters)
        save_sequence(sequence, output / f"sequence_n{n_vertices}.json")
        save_json(parameters, output / f"parameters_n{n_vertices}.json")
        save_json(
            {
                "n_vertices": n_vertices,
                "grid_positions": instance.grid_positions,
                "physical_coordinates_um": instance.coordinates(),
                "edges": instance.edges,
                "occupied_labels": list(sequence.register.qubit_ids),
                "layout_trap_count": len(sequence.register.layout.coords),
                "layout_coordinates_um": np.asarray(
                    sequence.register.layout.coords
                ).tolist(),
            },
            output / f"register_n{n_vertices}.json",
        )
        save_json(certificate, output / f"certificate_n{n_vertices}.json")
        report = {
            "n_vertices": n_vertices,
            "paper_target": instance.paper_target,
            "paper_valid_fraction": instance.paper_valid_fraction,
            "paper_500_shot_reproduction": {
                "replicates": paper_replicates,
                "mean_valid_approximation_ratio": float(
                    np.mean(
                        [
                            replicate["valid_approximation_ratio"]
                            for replicate in paper_replicates
                        ]
                    )
                ),
            },
            "baseline_exact": baseline_exact,
            "baseline_counts_500": baseline_500,
            "optimized_exact": optimized_exact,
            "optimized_counts_500": optimized_500,
            "optimized_modulated_exact": modulated_exact,
            "optimized_modulated_counts_500": modulated_500,
            "optimized_convergence": convergence,
            "robust_ensemble": robust_scores,
            "robust_worst_valid_approximation_ratio": min(
                entry["score"]["valid_approximation_ratio"] for entry in robust_scores
            ),
            "strictly_beats_paper": (
                modulated_500["valid_approximation_ratio"] > instance.paper_target
            ),
            "sequence": f"sequence_n{n_vertices}.json",
            "register": f"register_n{n_vertices}.json",
            "parameters": f"parameters_n{n_vertices}.json",
        }
        save_json(report, output / f"scores_n{n_vertices}.json")
        summary["instances"][str(n_vertices)] = report
    save_json(parameters, output / "parameters.json")
    summary["search_n11_score"] = optimized_n11
    save_json(summary, output / "scores.json")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Optimize exact Challenge 3 benchmarks.")
    parser.add_argument("--output", default="results/challenge03")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--hardware-retry",
        action="store_true",
        help="Optimize the 4 us N=17 schedule for a second hardware attempt.",
    )
    args = parser.parse_args()
    if args.hardware_retry:
        report = run_hardware_retry(args.output, quick=args.quick)
        print(
            "N=17 retry",
            f"modulated={report['modulated_exact']['valid_approximation_ratio']:.6f}",
            f"robust_worst={report['robust_worst_valid_approximation_ratio']:.6f}",
            f"ready={report['ready_for_hardware']}",
        )
        return
    report = run(args.output, quick=args.quick)
    for n_vertices, score in report["instances"].items():
        print(
            f"N={n_vertices}",
            f"paper={score['paper_target']:.6f}",
            "modulated_500="
            f"{score['optimized_modulated_counts_500']['valid_approximation_ratio']:.6f}",
            f"wins={score['strictly_beats_paper']}",
        )


if __name__ == "__main__":
    main()
