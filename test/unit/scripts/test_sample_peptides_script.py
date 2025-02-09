import random
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pyteomics")

from workflow.scripts.sample_peptides import (
    cleave_protein_sequence,
    convert_fasta_str_to_dict,
    get_accession_to_sequence_mapping,
    sample_accessions,
    sample_peptides,
    sample_peptides_from_sequence,
    UniProtConnector,
)


def test_get_fasta() -> None:
    fasta = UniProtConnector().get_fasta(["A0A3R6E0N5"])
    assert fasta == """\
>tr|A0A3R6E0N5|A0A3R6E0N5_9FIRM Glycosyl transferase OS=Roseburia intestinalis OX=166486 GN=DW264_18045 PE=4 SV=1
MCGGDDILKYRTYCKNQRDVAFVINGIIDEYWCGKLSEKEMKEDILTLYENNKEKLFKDG
QFTKIIQQQCGKKRINVISQILKNKLEKLE
"""


def test_get_fasta_returns_more_than_25_results() -> None:
    accessions = [
        'A0A3R6E0N5', 'A0A3E4LLW1', 'A0A7L6WN97', 'A0A6S6M774', 'A0A174UR68', 'A0A3R9NJY9', 'A0A0B0UEA8', 'A0A5C1JJV1',
        'A0A3E5GV22', 'A0A0L0LSS8', 'A0A427ZU75', 'A0A6H9PP43', 'A0A6G6L3K2', 'A0A6N3APR0', 'A0A0A1GRF5', 'A0A413EJW6',
        'A0A7I9AKF9', 'A0A1Q6B4J7', 'A0A558LVQ5', 'A0A2N5PYJ1', 'A0A139L2U9', 'F7LRI5', 'A0A376TL63', 'A0A8T3LEX4',
        'A0A139L7B9', 'A0A1S6GKC3', 'A0A7M1NVS6', 'A0A7D7DSC8', 'A0A0M1UJ95', 'A0A174PW89', 'A0A448PJ07', 'A0A415D293',
        'A0A7J4XQ64', 'A0A412P7C0', 'A0A0G3B8Y9', 'A0A9Q6F4E3', 'A0A7W9SEH0', 'A0A137SRT4', 'A0A2H4U0C6', 'A0A412RPZ4',
    ]
    fasta = UniProtConnector().get_fasta(accessions)
    assert fasta.count(">") > 25


def test_sample_peptides_from_sequence():
    expected = ['TGPNLHGLFGR', 'CSQCHTVEK', 'GIIWGEDTLME', 'TGQAPGYSYTAANK', 'IFIMKCSQCHTVEK']
    random.seed(12345)
    numpy_random_number_generator = np.random.default_rng(seed=random.randint(0, 2**31-1))
    sequence = "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLME"
    returned = sample_peptides_from_sequence(numpy_random_number_generator, sequence)
    assert returned == expected


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


def test_sample_peptides(workflow_path: Path) -> None:
    test_data = Path(__file__).parent.parent / "simulate_sample/sample_peptides/data/results/simulation/"
    sample_peptides(
        test_data / "tax2accessions.json",
        test_data / "sample_taxons_lineage_1.tsv",
        "/dev/null"
    )


def test_sample_accessions() -> None:
    tax_ids = ["1", "2", "3"]
    tax2acc = {
        "1": ["1a", "1b", "1c"],
        "2": ["2a", "2b", "2c"],
        "3": ["3a", "3b", "3c"],
    }
    random.seed(12345)
    returned = sample_accessions(tax2acc, tax_ids)
    assert returned == ['1b', '1a', '1c', '2b', '2a', '2c', '3c', '3b', '3a']


def test_get_accession_to_sequence_mapping():
    accessions = [
        'A0A3R6E0N5', 'A0A3E4LLW1', 'A0A7L6WN97', 'A0A6S6M774', 'A0A174UR68', 'A0A3R9NJY9', 'A0A0B0UEA8', 'A0A5C1JJV1',
        'A0A3E5GV22', 'A0A0L0LSS8', 'A0A427ZU75', 'A0A6H9PP43', 'A0A6G6L3K2', 'A0A6N3APR0', 'A0A0A1GRF5', 'A0A413EJW6',
        'A0A7I9AKF9', 'A0A1Q6B4J7', 'A0A558LVQ5', 'A0A2N5PYJ1', 'A0A139L2U9', 'F7LRI5', 'A0A376TL63', 'A0A8T3LEX4',
        'A0A139L7B9', 'A0A1S6GKC3', 'A0A7M1NVS6', 'A0A7D7DSC8', 'A0A0M1UJ95', 'A0A174PW89', 'A0A448PJ07', 'A0A415D293',
        'A0A7J4XQ64', 'A0A412P7C0', 'A0A0G3B8Y9', 'A0A9Q6F4E3', 'A0A7W9SEH0', 'A0A137SRT4', 'A0A2H4U0C6', 'A0A412RPZ4',
        'A0A8G1S9B2', 'A0A413I9L1', 'A0A2X6D6G0', 'O32560', 'A0A139KIL5', 'A0A9Q8DRL9', 'A0A250KGK4', 'A0A377DKM1',
        'A0A2N0UJJ5', 'A0A772E7H1', 'A0A3E5EFK6', 'A0A930LQ42', 'A0A827DYC8', 'A0A9P2HL15', 'A0A2N0URT8', 'A0A9Q4IVR9',
        'A0A4R6CTN2', 'A0A6P1Y4N3', 'A0A8B3B8D5', 'A0A3R5ZQX8', 'A0A174ERA0', 'A0A6D0FWW0', 'J1LM39', 'A0A7M1NXV3',
        'A0A6A1Z4W4', 'A0A3S4K6R1', 'A0A1V3CDT4', 'A0A4Q5HEZ6', 'A0A8T5ZQ29', 'A0A3E5ET67', 'A0A4Q5IDM6', 'A0A413YX76',
        'A0A1V2FXF9', 'A0A5M6A476', 'A0A2U2RUX8', 'A0A849YN47', 'A0A2A3UQ42', 'A0A0L0LST0', 'A0A828ABF2', 'A0A826J2R3',
        'A0A376K0H5', 'A0A174NBZ0', 'A0A139LVH8', 'A0A930Q174', 'A0A2N0SPF6', 'A0A1Q8I2P6', 'A0A2X1JGP1', 'A0A2X3K270',
        'A0A7J4XUU1', 'A0A173RK00', 'A0A827QNQ9', 'A0A9Q7ZIA1', 'A0A378TZE8', 'A0A921K5C8'
    ]
    uniprot_connector = UniProtConnector()
    assert len(get_accession_to_sequence_mapping(accessions, uniprot_connector)) > 25
