import random
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pyteomics")

from workflow.scripts.sample_peptides import sample_peptides_from_sequence, cleave_protein_sequence, \
    convert_fasta_str_to_dict, get_protein_sequence, sample_peptides


def test_sample_peptides_from_sequence():
    random.seed(12345)
    numpy_random_number_generator = np.random.default_rng(seed=random.getstate()[1][0])
    sequence = "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLME"
    assert (sample_peptides_from_sequence(numpy_random_number_generator, sequence) ==
            ['CSQCHTVEK', 'TGPNLHGLFGR', 'GIIWGEDTLME', 'TGQAPGYSYTAANK', 'IFIMKCSQCHTVEK'])


@pytest.mark.parametrize(
    ("missed_cleavages", "expected"),
    [
        (0, {'TGQAPGYSYTAANK', 'TGPNLHGLFGR', 'GIIWGEDTLME', 'CSQCHTVEK'}),
        (1, {'CSQCHTVEKGGK', 'KTGQAPGYSYTAANK', 'HKTGPNLHGLFGR', 'TGQAPGYSYTAANKNK', 'TGPNLHGLFGRK', 'IFIMKCSQCHTVEK', 'NKGIIWGEDTLME', 'MGDVEKGK'}),
        (2, {'KTGQAPGYSYTAANKNK', 'GKKIFIMK', 'TGPNLHGLFGRKTGQAPGYSYTAANK', 'TGQAPGYSYTAANKNKGIIWGEDTLME', 'KIFIMKCSQCHTVEK', 'IFIMKCSQCHTVEKGGK', 'GGKHKTGPNLHGLFGR', 'HKTGPNLHGLFGRK', 'MGDVEKGKK', 'CSQCHTVEKGGKHK'}),
    ]
)
def test_cleave_protein_sequence(missed_cleavages: int, expected: set[str]):
    sequence = "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLME"
    assert cleave_protein_sequence(sequence, missed_cleavages=missed_cleavages) == expected


def test_fasta_to_dict():
    fasta = """>sp|P12345|AATM_RABIT Aspartate aminotransferase, mitochondrial OS=Oryctolagus cuniculus OX=9986 GN=GOT2 PE=1 SV=2
MALLHSARVLSGVASAFHPGLAAAASARASSWWAHVEMGPPDPILGVTEAYK
>sp|P99999|CYC_HUMAN Cytochrome c OS=Homo sapiens OX=9606 GN=CYCS PE=1 SV=2
MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIW
GEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE"""
    assert convert_fasta_str_to_dict(fasta) == {
        "P12345": "MALLHSARVLSGVASAFHPGLAAAASARASSWWAHVEMGPPDPILGVTEAYK",
        "P99999": "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE",
    }


def test_get_protein_sequence():
    assert (
        get_protein_sequence("P99999")
        == "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE"
    )


def test_sample_peptides(workflow_path: Path) -> None:
    sample_peptides(
        Path(__file__).parent.parent / "simulate_sample/sample_peptides/data/results/simulation/tax2accessions_1.json",
        "/dev/null"
    )
