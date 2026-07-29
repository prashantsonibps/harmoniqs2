from __future__ import annotations

import numpy as np

from harmoniqs.challenge01 import (
    REFERENCE_DURATION_NS,
    REFERENCE_OMEGA,
    bell_fidelity_symmetric,
    evolve_controls,
    reference_sequence,
    simulate_sequence,
)
from harmoniqs.challenge02 import (
    GRAPH_EDGES,
    REFERENCE,
    make_sequence,
    register_coordinates,
)


def pairwise_distances(coordinates: list[tuple[float, float]]) -> np.ndarray:
    coords = np.asarray(coordinates)
    return np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)


def test_challenge01_dense_model_matches_pulser() -> None:
    spacing = 5.0
    model = evolve_controls(
        np.array([[REFERENCE_OMEGA, 0.0]]),
        spacing,
        REFERENCE_DURATION_NS,
    )
    model_fidelity = bell_fidelity_symmetric(model)
    pulser_fidelity, _ = simulate_sequence(reference_sequence(spacing))
    assert abs(model_fidelity - pulser_fidelity) < 2e-4


def test_star_embedding_has_only_center_leaf_edges() -> None:
    coords = register_coordinates("star", 5.5)
    distances = pairwise_distances(coords)
    radius = 7.2
    found = {
        (i, j)
        for i in range(4)
        for j in range(i + 1, 4)
        if distances[i, j] < radius
    }
    assert found == set(GRAPH_EDGES["star"])


def test_c5_embedding_has_only_cycle_edges() -> None:
    coords = register_coordinates("c5", 5.5)
    distances = pairwise_distances(coords)
    radius = 7.2
    found = {
        (i, j)
        for i in range(5)
        for j in range(i + 1, 5)
        if distances[i, j] < radius
    }
    assert found == set(GRAPH_EDGES["c5"])


def test_challenge02_sequences_fit_device() -> None:
    for graph in GRAPH_EDGES:
        sequence = make_sequence(graph, REFERENCE)
        assert sequence.get_duration() == 4000
