from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Iterable, Sequence as TypingSequence

import numpy as np
from pulser import AnalogDevice, Sequence

SHOTS = 500
C6_OVER_HBAR = float(AnalogDevice.interaction_coeff)
OMEGA_MAX = float(AnalogDevice.channels["rydberg_global"].max_amp)
DETUNING_MAX = float(AnalogDevice.channels["rydberg_global"].max_abs_detuning)
CLOCK_NS = int(AnalogDevice.channels["rydberg_global"].clock_period)
MAX_DURATION_NS = int(AnalogDevice.max_sequence_duration)


def device_summary() -> dict[str, float | int]:
    channel = AnalogDevice.channels["rydberg_global"]
    return {
        "interaction_coeff_rad_us_um6": float(AnalogDevice.interaction_coeff),
        "max_atoms": int(AnalogDevice.max_atom_num),
        "min_atom_distance_um": float(AnalogDevice.min_atom_distance),
        "max_radius_um": float(AnalogDevice.max_radial_distance),
        "max_duration_ns": int(AnalogDevice.max_sequence_duration),
        "clock_period_ns": int(channel.clock_period),
        "max_omega_rad_us": float(channel.max_amp),
        "max_abs_detuning_rad_us": float(channel.max_abs_detuning),
        "max_runs": int(AnalogDevice.max_runs),
    }


def validate_sequence(sequence: Sequence) -> None:
    sequence.get_duration()
    AnalogDevice.validate_register(sequence.register)
    if sequence.get_duration() > MAX_DURATION_NS:
        raise ValueError("Sequence exceeds the AnalogDevice duration limit.")


def save_sequence(sequence: Sequence, path: str | Path) -> Path:
    validate_sequence(sequence)
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(sequence.to_abstract_repr(), encoding="utf-8")
    return destination


def save_json(payload: object, path: str | Path) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return destination


def probabilities_from_state(
    state: np.ndarray, *, pulser_basis: bool = False
) -> dict[str, float]:
    """Return probabilities with ``1`` consistently meaning Rydberg-excited.

    QuTiP/Pulser orders each local basis as ``|r>, |g>`` while the dense
    optimizers use the conventional binary order ``|g>, |r>``.
    """
    vector = np.asarray(state, dtype=complex).reshape(-1)
    n_qubits = int(np.log2(vector.size))
    probabilities = {}
    for index, amplitude in enumerate(vector):
        bitstring = format(index, f"0{n_qubits}b")
        if pulser_basis:
            bitstring = "".join("1" if bit == "0" else "0" for bit in bitstring)
        probabilities[bitstring] = float(abs(amplitude) ** 2)
    return probabilities


def deterministic_counts(
    probabilities: dict[str, float], shots: int = SHOTS
) -> Counter[str]:
    """Convert probabilities to exactly ``shots`` counts without RNG noise."""
    expected = {key: shots * value for key, value in probabilities.items()}
    counts = Counter({key: int(value) for key, value in expected.items()})
    remainder = shots - sum(counts.values())
    order = sorted(
        expected,
        key=lambda key: expected[key] - counts[key],
        reverse=True,
    )
    counts.update(order[:remainder])
    return +counts


def is_independent(bitstring: str, edges: Iterable[tuple[int, int]]) -> bool:
    return all(not (bitstring[i] == bitstring[j] == "1") for i, j in edges)


def mis_bitstrings(
    n_vertices: int, edges: TypingSequence[tuple[int, int]]
) -> set[str]:
    independent = [
        format(value, f"0{n_vertices}b")
        for value in range(2**n_vertices)
        if is_independent(format(value, f"0{n_vertices}b"), edges)
    ]
    maximum = max(bits.count("1") for bits in independent)
    return {bits for bits in independent if bits.count("1") == maximum}


def mis_probability(
    probabilities: dict[str, float], edges: TypingSequence[tuple[int, int]]
) -> float:
    targets = mis_bitstrings(len(next(iter(probabilities))), edges)
    return float(sum(probabilities.get(bits, 0.0) for bits in targets))
