from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pulser import Sequence

from harmoniqs.challenge01 import (
    REFERENCE_DURATION_NS,
    REFERENCE_OMEGA,
    bell_fidelity_symmetric,
    evolve_controls,
    reference_sequence,
    simulate_sequence,
)

SHAPED_DIR = Path(__file__).resolve().parents[1] / "results" / "challenge01" / "shaped"

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


def test_shaped_pulse_scores_exist() -> None:
    scores_path = SHAPED_DIR / "scores.json"
    assert scores_path.is_file(), f"Missing {scores_path}"
    data = json.loads(scores_path.read_text(encoding="utf-8"))
    assert data["challenge"] == 1
    for key in ("5.0um", "6.5um"):
        spacing = data["spacings"][key]
        f_pulser_4ns = spacing["optimized_fidelity_pulser_4ns"]
        assert 0.99 < f_pulser_4ns < 1.0
        assert spacing["within_device_bounds"]
        assert spacing["strictly_beats_baseline"]


def test_shaped_pulse_sequences_are_valid_pulser() -> None:
    for key, spacing in [("5.0um", 5.0), ("6.5um", 6.5)]:
        seq_path = SHAPED_DIR / f"sequence_{key}.json"
        assert seq_path.is_file(), f"Missing {seq_path}"
        representation = seq_path.read_text(encoding="utf-8")
        sequence = Sequence.from_abstract_repr(representation)
        duration = sequence.get_duration()
        assert duration % 4 == 0, f"{key} duration {duration} not divisible by 4 ns clock"
        assert duration <= 6000, f"{key} duration {duration} exceeds device limit"


def test_shaped_pulse_cross_validation_consistent() -> None:
    cv_path = SHAPED_DIR / "cross_validation_results.json"
    assert cv_path.is_file(), f"Missing {cv_path}"
    data = json.loads(cv_path.read_text(encoding="utf-8"))
    for run in ("r1_5um", "r2_6.5um"):
        entry = data[run]
        qt_4ns = entry["F_qt_4ns"]
        pulser_4ns = entry["F_pulser_4ns"]
        diff = abs(qt_4ns - pulser_4ns)
        assert diff < 2e-7, (
            f"{run}: QuTiP ({qt_4ns:.10f}) vs Pulser ({pulser_4ns:.10f}) "
            f"disagree by {diff:.2e}"
        )
        assert entry["constraints"]["amp_ok"]
        assert entry["constraints"]["det_ok"]


def test_shaped_pulse_submission_exists() -> None:
    sub_path = SHAPED_DIR / "submission.json"
    assert sub_path.is_file(), f"Missing {sub_path}"
    data = json.loads(sub_path.read_text(encoding="utf-8"))
    assert data["challenge"] == 1
    assert data["team"] == "team 2"
    assert "explanation" in data
    assert len(data["explanation"]) > 50
    for spacing_key in ("r1_5um", "r2_6.5um"):
        summary = data["summary"][spacing_key]
        assert summary["optimized_fidelity"] > 0.999
        assert summary["improvement"] > 0
