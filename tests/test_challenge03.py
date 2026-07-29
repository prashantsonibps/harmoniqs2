from __future__ import annotations

import numpy as np
from pulser import Sequence

from harmoniqs.challenge02 import _operators
from harmoniqs.challenge03 import (
    SparseRydbergModel,
    certify_instance,
    evolve_sparse,
    load_instance,
    make_sequence,
    modulated_controls,
    PAPER_500_SHOT_REPLICATES,
    PAPER_TARGETS,
    paper_reference_parameters,
    probabilities_from_state,
    qubo_cost,
    score_counts,
)


def test_exact_paper_instances_and_certificates() -> None:
    expected = {11: (17, 4, 21), 13: (27, 4, 30), 17: (31, 6, 24)}
    for n_vertices, (edge_count, alpha, optimum_count) in expected.items():
        instance = load_instance(n_vertices)
        certificate = certify_instance(instance)
        assert len(instance.edges) == edge_count
        assert certificate["alpha"] == alpha
        assert certificate["minimum_qubo_cost"] == -alpha
        assert certificate["number_of_maximum_independent_sets"] == optimum_count


def test_paper_scoring_is_conditioned_on_valid_samples() -> None:
    instance = load_instance(11)
    counts = {
        "10011100000": 2,  # valid optimum, size 4
        "10000000000": 1,  # valid size 1
        "11000000000": 1,  # invalid
    }
    score = score_counts(counts, instance)
    assert score["shots"] == 4
    assert np.isclose(score["valid_probability"], 0.75)
    assert np.isclose(score["valid_approximation_ratio"], 9 / 12)
    assert qubo_cost("11000000000", instance.edges) == 0


def test_published_scores_reconstruct_from_500_shot_replicates() -> None:
    for n_vertices, replicates in PAPER_500_SHOT_REPLICATES.items():
        alpha = certify_instance(load_instance(n_vertices))["alpha"]
        ratios = [
            replicate["valid_size_sum"] / (alpha * replicate["valid_count"])
            for replicate in replicates
        ]
        assert np.isclose(np.mean(ratios), PAPER_TARGETS[n_vertices], atol=5e-8)


def test_sparse_model_matches_dense_operators() -> None:
    coordinates = [(0.0, 0.0), (5.0, 0.0), (0.0, 6.0)]
    sparse = SparseRydbergModel.from_coordinates(coordinates)
    drive, excitations, interactions = _operators(coordinates)
    assert np.allclose(sparse.drive.toarray(), drive)
    assert np.allclose(sparse.excitations, excitations)
    assert np.allclose(sparse.interactions, interactions)


def test_sparse_state_is_normalized() -> None:
    instance = load_instance(11)
    model = SparseRydbergModel.from_coordinates(instance.coordinates())
    parameters = paper_reference_parameters()
    state = evolve_sparse(model, parameters, steps=8)
    probabilities = probabilities_from_state(state, instance.n_vertices)
    assert np.isclose(sum(probabilities.values()), 1.0, atol=1e-8)


def test_sequences_have_sixty_trap_layouts_and_round_trip() -> None:
    for n_vertices in (11, 13, 17):
        sequence = make_sequence(load_instance(n_vertices), paper_reference_parameters())
        assert len(sequence.register.layout.coords) == 60
        assert sequence.get_duration() == 4000
        restored = Sequence.from_abstract_repr(sequence.to_abstract_repr())
        assert len(restored.register.qubit_ids) == n_vertices
        assert len(restored.register.layout.coords) == 60


def test_modulation_sampling_includes_channel_tail() -> None:
    sequence = make_sequence(load_instance(11), paper_reference_parameters())
    amplitude, detuning, durations = modulated_controls(sequence)
    assert len(amplitude) == len(detuning) == len(durations)
    assert np.isclose(durations.sum(), 4.168)
