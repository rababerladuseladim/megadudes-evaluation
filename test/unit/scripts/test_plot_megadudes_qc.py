from workflow.scripts.plot_megadudes_qc import read_ground_truth_file, get_diamond_hit_counts, get_value_overlap
from pathlib import Path
from pandas.testing import assert_frame_equal, assert_series_equal
import pandas as pd
import pytest

TEST_DATA = Path(__file__).parent / "plot_megadudes_qc"


def test_read_ground_truth_file() -> None:
    expected = pd.DataFrame(
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
    returned = read_ground_truth_file(TEST_DATA / "ground_truth.csv")
    assert_frame_equal(expected, returned)


# def test_get_diamond_hit_counts() -> None:
#     expected = pd.DataFrame()
#     returned = get_diamond_hit_counts()

@pytest.mark.parametrize(
    ("s1", "s2", "expected"),
    (
        (
            pd.Series([1,2,3]),
            pd.Series([2,3,4]),
            pd.Series({1: "FP", 2: "TP", 3: "TP", 4: "FN"})
        ),
        (
            pd.Series([1,2,3,pd.NA]),
            pd.Series([2,3,4,"-"]),
            pd.Series({1: "FP", 2: "TP", 3: "TP", 4: "FN"})
        ),
    )
)
def test_get_value_overlap(s1, s2, expected) -> None:
    assert_series_equal(expected, get_value_overlap(s1, s2))

# def test_get_diamond_hit_counts
# def test_get_unipept_hit_counts
# def test_get_megadudes_hit_counts
# def test_calc_eval_metrics
