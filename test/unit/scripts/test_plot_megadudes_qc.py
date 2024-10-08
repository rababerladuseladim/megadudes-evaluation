from workflow.scripts.plot_megadudes_qc import (
    calc_eval_metrics,
    read_ground_truth_file,
    get_diamond_hit_counts,
    get_value_overlap,
    get_unipept_hit_counts,
    get_megadudes_hit_counts,
)
from pathlib import Path
from pandas.testing import assert_frame_equal, assert_series_equal
import pandas as pd
import numpy as np
import pytest

TEST_DATA = Path(__file__).parent / "plot_megadudes_qc"


@pytest.fixture()
def ground_truth_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "superkingdom": {0: 2, 1: 2},
            "phylum": {0: 1224, 1: 1224},
            "class": {0: 28211, 1: 1236},
            "order": {0: 356, 1: 135622},
            "family": {0: 82115, 1: 72275},
            "genus": {0: 357, 1: 226},
            "species": {0: 358, 1: 28108},
            "subspecies": {0: pd.NA, 1: 1004788},
        },
        dtype=pd.Int64Dtype(),
    )


def test_read_ground_truth_file(ground_truth_df: pd.DataFrame) -> None:
    returned = read_ground_truth_file(TEST_DATA / "ground_truth.csv")
    assert_frame_equal(ground_truth_df, returned)


@pytest.mark.parametrize(
    ("s1", "s2", "expected"),
    (
        (
            pd.Series([1, 2, 3]),
            pd.Series([2, 3, 4]),
            pd.Series({1: "FP", 2: "TP", 3: "TP", 4: "FN"}),
        ),
        (
            pd.Series([1, 2, 3, pd.NA]),
            pd.Series([2, 3, 4, "-"]),
            pd.Series({1: "FP", 2: "TP", 3: "TP", 4: "FN"}),
        ),
    ),
)
def test_get_value_overlap(s1, s2, expected) -> None:
    assert_series_equal(expected, get_value_overlap(s1, s2))


def test_get_diamond_hit_counts(ground_truth_df: pd.DataFrame) -> None:
    expected = pd.DataFrame(
        {
            "superkingdom": {"TP": 1.0, "FP": 0.0, "FN": 0.0},
            "phylum": {"TP": 1.0, "FP": 3.0, "FN": 0.0},
            "class": {"TP": 1, "FP": 2, "FN": 1},
            "order": {"TP": 0.0, "FP": 4.0, "FN": 2.0},
            "family": {"TP": 0.0, "FP": 5.0, "FN": 2.0},
            "genus": {"TP": 0.0, "FP": 5.0, "FN": 2.0},
            "species": {"TP": 0.0, "FP": 9.0, "FN": 2.0},
            "subspecies": {"TP": 0.0, "FP": 0.0, "FN": 1.0},
            "eval": {"TP": "TP", "FP": "FP", "FN": "FN"},
        }
    )
    returned = get_diamond_hit_counts(TEST_DATA / "diamond.tsv", ground_truth_df)
    assert_frame_equal(expected, returned)


def test_get_unipept_hit_counts(ground_truth_df: pd.DataFrame) -> None:
    expected = pd.DataFrame(
        {
            "superkingdom": {"TP": 1.0, "FP": 0.0, "FN": 0.0},
            "phylum": {"TP": 0.0, "FP": 1.0, "FN": 1.0},
            "class": {"TP": 0.0, "FP": 1.0, "FN": 2.0},
            "order": {"TP": 0.0, "FP": 1.0, "FN": 2.0},
            "family": {"TP": 0.0, "FP": 1.0, "FN": 2.0},
            "genus": {"TP": 0.0, "FP": 1.0, "FN": 2.0},
            "species": {"TP": 0.0, "FP": 1.0, "FN": 2.0},
            "subspecies": {"TP": 0.0, "FP": 0.0, "FN": 1.0},
            "eval": {"TP": "TP", "FP": "FP", "FN": "FN"},
        }
    )
    returned = get_unipept_hit_counts(TEST_DATA / "unipept.tsv", ground_truth_df)
    assert_frame_equal(expected, returned)


def test_get_megadudes_hit_counts(ground_truth_df: pd.DataFrame) -> None:
    expected = pd.DataFrame(
        {
            "superkingdom": {"TP": 1, "FP": 0, "FN": 0},
            "phylum": {"TP": 0, "FP": 3, "FN": 1},
            "class": {"TP": 0, "FP": 0, "FN": 2},
            "order": {"TP": 0, "FP": 0, "FN": 2},
            "family": {"TP": 0, "FP": 0, "FN": 2},
            "genus": {"TP": 0, "FP": 0, "FN": 2},
            "species": {"TP": 0, "FP": 0, "FN": 2},
            "subspecies": {"TP": 0, "FP": 0, "FN": 1},
            "eval": {"TP": "TP", "FP": "FP", "FN": "FN"},
        }
    )
    returned = get_megadudes_hit_counts(TEST_DATA / "megadudes.tsv", ground_truth_df)
    assert_frame_equal(expected, returned)


@pytest.mark.parametrize(
    ("df_hits", "expected"),
    [
        (
            pd.DataFrame(
                {
                    "superkingdom": {"TP": 1, "FP": 1, "FN": 1},
                    "method": {"TP": "megadudes", "FP": "megadudes", "FN": "megadudes"},
                    "eval": {"TP": "TP", "FP": "FP", "FN": "FN"},
                }
            ),
            pd.DataFrame(
                {
                    "method": {0: "megadudes"},
                    "taxon_level": {0: "superkingdom"},
                    "sensitivity": {0: 50.0},
                    "precision": {0: 50.0},
                    "F1-score": {0: 50.0},
                    "fdr": {0: 50.0},
                }
            ),
        ),
        (
            pd.DataFrame(
                {
                    "superkingdom": {0: 1, 1: 1, 2: 1, 3: 3, 4: 1, 5: 1},
                    "method": {0: "megadudes", 1: "megadudes", 2: "megadudes", 3: "diamond", 4: "diamond", 5: "diamond"},
                    "eval": {0: "TP", 1: "FP", 2: "FN", 3: "TP", 4: "FP", 5: "FN"},
                }
            ),
            pd.DataFrame(
                {
                    "method": {0: "megadudes", 1: "diamond"},
                    "taxon_level": {0: "superkingdom", 1: "superkingdom"},
                    "sensitivity": {0: 50.0, 1: 75.0},
                    "precision": {0: 50.0, 1: 75.0},
                    "F1-score": {0: 50.0, 1: 75.0},
                    "fdr": {0: 50.0, 1: 25.0},
                }
            ),
        ),
        (
            pd.DataFrame(
                {
                    "superkingdom": {"TP": 0, "FP": 0, "FN": 1},
                    "method": {"TP": "megadudes", "FP": "megadudes", "FN": "megadudes"},
                    "eval": {"TP": "TP", "FP": "FP", "FN": "FN"},
                }
            ),
            pd.DataFrame(
                {
                    "method": {0: "megadudes"},
                    "taxon_level": {0: "superkingdom"},
                    "sensitivity": {0: 0.0},
                    "precision": {0: np.nan},
                    "F1-score": {0: np.nan},
                    "fdr": {0: np.nan},
                }
            ),
        ),
    ],
    ids=[
        "single_method",
        "multi_method",
        "no_positives"
    ]
)
def test_calc_eval_metrics(df_hits: pd.DataFrame, expected: pd.DataFrame) -> None:
    returned = calc_eval_metrics(df_hits)
    assert_frame_equal(expected, returned)
