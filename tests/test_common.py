from __future__ import annotations

import numpy as np
from pulser import Sequence

from harmoniqs.challenge01 import reference_sequence
from harmoniqs.common import (
    deterministic_counts,
    device_summary,
    mis_bitstrings,
    mis_probability,
)


def test_device_envelope_matches_challenge() -> None:
    summary = device_summary()
    assert summary["max_atoms"] == 80
    assert summary["min_atom_distance_um"] == 5.0
    assert summary["max_duration_ns"] == 6000
    assert summary["clock_period_ns"] == 4


def test_mis_scoring() -> None:
    star_edges = ((0, 1), (0, 2), (0, 3))
    assert mis_bitstrings(4, star_edges) == {"0111"}
    probabilities = {"0111": 0.7, "1000": 0.2, "0000": 0.1}
    assert np.isclose(mis_probability(probabilities, star_edges), 0.7)


def test_deterministic_counts_preserves_shots() -> None:
    counts = deterministic_counts({"00": 0.333, "01": 0.333, "10": 0.334})
    assert sum(counts.values()) == 500


def test_sequence_round_trip() -> None:
    sequence = reference_sequence(5.0)
    restored = Sequence.from_abstract_repr(sequence.to_abstract_repr())
    assert restored.get_duration() == 352
    assert len(restored.register.qubit_ids) == 2
