from workflow.scripts.plot_megadudes_qc import read_ground_truth_file, get_diamond_hit_counts, get_value_overlap, \
    TAX_LEVELS
from pathlib import Path
from pandas.testing import assert_frame_equal, assert_series_equal
import pandas as pd
import pytest

TEST_DATA = Path(__file__).parent / "plot_megadudes_qc"


@pytest.fixture()
def ground_truth_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            'superkingdom': {0: 2, 1: 2},
            'phylum': {0: 1224, 1: 1224},
            'class': {0: 28211, 1: 1236},
            'order': {0: 356, 1: 135622},
            'family': {0: 82115, 1: 72275},
            'genus': {0: 357, 1: 226},
            'species': {0: 358, 1: 28108},
            'subspecies': {0: pd.NA, 1: 1004788}
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
            pd.Series({1: "FP", 2: "TP", 3: "TP", 4: "FN"})
        ),
        (
            pd.Series([1, 2, 3, pd.NA]),
            pd.Series([2, 3, 4, "-"]),
            pd.Series({1: "FP", 2: "TP", 3: "TP", 4: "FN"})
        ),
    )
)
def test_get_value_overlap(s1, s2, expected) -> None:
    assert_series_equal(expected, get_value_overlap(s1, s2))


def test_get_diamond_hit_counts(ground_truth_df: pd.DataFrame) -> None:
    expected = pd.DataFrame(
        {
            'superkingdom': {'TP': 1.0, 'FP': 0.0, 'FN': 0.0},
            'phylum': {'TP': 1.0, 'FP': 3.0, 'FN': 0.0},
            'class': {'TP': 1, 'FP': 2, 'FN': 1},
            'order': {'TP': 0.0, 'FP': 4.0, 'FN': 2.0},
            'family': {'TP': 0.0, 'FP': 5.0, 'FN': 2.0},
            'genus': {'TP': 0.0, 'FP': 5.0, 'FN': 2.0},
            'species': {'TP': 0.0, 'FP': 9.0, 'FN': 2.0},
            'subspecies': {'TP': 0.0, 'FP': 0.0, 'FN': 1.0},
            'eval': {'TP': 'TP', 'FP': 'FP', 'FN': 'FN'},
            'method': {'TP': 'diamond', 'FP': 'diamond', 'FN': 'diamond'}
        }
    )
    returned = get_diamond_hit_counts(TEST_DATA / "diamond.tsv", ground_truth_df)
    assert_frame_equal(expected, returned)

# def test_get_unipept_hit_counts
# def test_get_megadudes_hit_counts
# def test_calc_eval_metrics
