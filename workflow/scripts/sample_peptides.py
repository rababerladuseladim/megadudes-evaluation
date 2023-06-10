import random
from itertools import groupby
from pathlib import Path
from typing import cast, Dict, Iterable
from pyteomics.parser import cleave
import requests


import json
import sys

LOG_HANDLE = sys.stderr


def sample_peptides(accessions_file, output):
    with open(accessions_file, "r") as handle:
        tax2acc = cast(Dict[str, list[str]], json.load(handle))

    # sample accessions
    random.seed(str(output))
    accessions_sample = []
    for i, (tax_id, accessions) in enumerate(tax2acc.items()):
        size = min(len(accessions), 100)
        accessions_sample.extend(random.sample(accessions, size))
        if i >= 100:
            break

    print(f"Sampled {len(accessions_sample)} accessions", file=LOG_HANDLE)
    seed = random.getstate()

    # get protein sequences for sampled accessions
    fastas = []
    chunk_size = 100
    connector = UniProtConnector()
    for chunk_start in range(0, len(accessions_sample), chunk_size):
        chunk = accessions_sample[chunk_start:chunk_start + chunk_size]
        fastas.append(connector.get_fasta(chunk))
    acc2seq = {acc: seq for f in fastas for acc, seq in convert_fasta_str_to_dict(f).items()}

    # ensure deterministic behaviour although previous fetching of sequences takes a non-deterministic amount of
    # computation
    random.setstate(seed)

    # cleave proteins sequences and sample peptides
    sequences = sorted(acc2seq.values())
    peptides_sample = list()
    for seq in sequences:
        peptides = sorted(cleave_protein_sequence(seq))
        peptides_sample.extend(random.sample(list(peptides), min(3, len(peptides))))

    print(f"Sampled {len(peptides_sample)} peptides", file=LOG_HANDLE)

    # write fasta
    with open(output, "w") as handle:
        handle.write("\n".join(peptides_sample) + "\n")


def cleave_protein_sequence(sequence: str, min_length: int = 5):
    return cleave(sequence, rule="trypsin", min_length=min_length)


def get_protein_sequence(accession):
    connector = UniProtConnector()
    fasta = connector.get_fasta([accession])
    return "".join(fasta.split("\n")[1:])


class UniProtConnector:
    url = "https://rest.uniprot.org/uniprotkb/"

    def __init__(self):
        self.session = requests.Session()
        self.uniprot_version = requests.get(self.url).headers.get('X-UniProt-Release')
        print(f"Uniprot Release Number: {self.uniprot_version}", file=LOG_HANDLE)

    def get_fasta(self, uniprot_accessions: Iterable[str]):
        uniprot_accessions = ["accession%3A" + acc for acc in uniprot_accessions]
        url = self.url + f"search?query={'+OR+'.join(uniprot_accessions)}&format=fasta"
        response = self.session.get(url)
        response.raise_for_status()
        fasta = response.text
        return fasta


def convert_fasta_str_to_dict(fasta):
    protein_dict = {}
    fasta = fasta.strip()
    fasta_iter = (x[1] for x in groupby(fasta.split("\n"), lambda line: line[0] == ">"))

    for header in fasta_iter:
        # keep only uniprot accession, between pipe-chars ('|')
        acc = header.__next__().split("|")[1]

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in fasta_iter.__next__())

        protein_dict[acc] = seq
    return protein_dict


def test_cleave_protein_sequence():
    seq = "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE"
    peps = {
        "GIIWGEDTLMEYLENPK",
        "TGQAPGYSYTAANK",
        "TGPNLHGLFGR",
        "CSQCHTVEK",
        "ADLIAYLK",
        "MIFVGIK",
        "YIPGTK",
        "MGDVEK",
        "IFIMK",
    }
    assert cleave_protein_sequence(seq) == peps


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


def test_sample_peptides(tmpdir):
    sample_peptides(
        Path(__file__).parent.parent.parent / "test/unit/simulate_sample/sample_peptides/data/results/map_taxids_to_uniprot_accessions/tax2accessions.json",
        "/dev/null"
        # tmpdir / "peptides.txt"
    )
    expected_accestions = "P68254 Q9CQV8 P63101 P68510 P61982 P62259 O70456 P62260 P63102 P68511 P35213 P68255 P61983 Q9MB95 Q80HU0 P18558 P18559 P18560 P68744 Q65209 Q65210 P18557 P0C9H2 P0C9K1 P0C9F1 P0C9J9 P0C9I7 P0C9H5 P0C9F9 P0C9J5 P0C9I3 P0C9G7 P0C9F5 P0C9K3 P0C9G3 P0C9H9 P16536 P05080 P21215 Q8GBW6"
    expected_peptides = "IASALEMYATK LNNLLDPHSFDEVGAFR IVDWGDYLEVK QIFLTAIGDQAR GDSVHALCALWK GFPCK EAAENSLVAYK GDYHR YDEMVESMK DSTLIMQLLR IISSIEQK GDYHR SVTEQGAELSNEER DSTLIMQLLR GIVDQSQQAYQEAFEISK TEGAEK QQMAR DNLTLWTSDTQGDEAEAGEGGEN QTIENSQGAYQEAFDISK LGLALNFSVFYYEILNNPELACTLAK YLAEVACGDDR YDDMATCMK NVVGGR YLAEVACGDDR YLAEVATGDDK SAYQEAMDISK VETELR NLLSVAYK AVTELNEPLSNEDR QAFDDAIAELDTLNEDSYK QAFDDAIAELDTLNEDSYK NLLSVAYK GDYYR QYCLYFIIGIAYTDCFICALCK MGGGGDHQQLSIK VLDDMPLIVQNDYISK NVAIFQDYHGFPEFR CVEPGWFR GLEEVGINCLK CMYEAHFR DHPIITQNDYIVNCTVSR LLALLSILIWLSQPALNRPLSIFYMK CSYEK ELEHWCTHGK NHSPILENNYIANCSIYR LLTYGFYLVGCVLVANYVR IGSIPPER SLQTIMYLMVLLVIFFLLSQLMLYR ITIDVTPK DVLLAQSVAVEEAK IGNVVTISYNLEK MLVIFLGILGLLANQVLGLPIQAGGHLCSTDNPPQEELGYWCTYMESCK VNESMPLIIENSYLTSCEVSR WYNQCTYGEGNGHYHVMDCSNPVPHNRPHQLR WYNQCTYSEGNGHYHVMDCSNPVPHNRPHR FCWECAHGICK MLVIFLGILGLLANQVLGLPTQAGGHLR MLVIFLGILGLLANQVSSQLVGQLHPTENPSENELEYWCTYMECCQFCWDCQNGLCVNK LGNTTILENEYVHPCIVSR CMYDLDK IDGSAIYK NEYVK QTTVSNSQQAYQEAFEISK ACSLAK DSTLIMQLLR VISSIEQK IEAELQDICSDVLELLDK LAEQAER ELEAVCQDVLSLLDNYLIK VFYLK LAEQAER YDDMAAAMK NVVGAR TSADGNEK"


if "snakemake" in globals():
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        sample_peptides(
            accessions_file=snakemake.input["accessions"], output=snakemake.output[0]
        )
