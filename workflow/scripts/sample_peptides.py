from itertools import groupby
from pyteomics.parser import cleave
from urllib import request


import sys

LOG_HANDLE = sys.stderr


def sample_peptides(accessions_file, output):
    pass


def cleave_protein_sequence(sequence, min_length=5):
    return cleave(sequence, rule="trypsin", min_length=min_length)


def get_protein_sequence(accession):
    fasta = get_protein_sequences([accession])
    return "".join(fasta.split("\n")[1:])


def get_protein_sequences(uniprot_list):
    """Retrieves the sequences from the UniProt database based on the list of
    UniProt ids.
    In general,
        1. Compose your query here with the advanced search tool:
    https://www.uniprot.org/uniprot/?query=id%3Ap40925+OR+id%3Ap40926+OR+id%3Ao43175&sort=score
        2. Replace `&sort=score` with `&format=fasta`
        3. Edit this function as necessary
    Returns:
        protein_dict (dict): the updated dictionary
    """
    # This makes it so we match only the ENTRY field
    uniprot_list = ["id%3A" + id for id in uniprot_list]
    line = "+OR+".join(uniprot_list)
    url = f"https://www.uniprot.org/uniprot/?query={line}&format=fasta"
    with request.urlopen(url) as f:
        fasta = f.read().decode("utf-8").strip()
    return fasta


def fasta_str_to_dict(fasta):
    protein_dict = {}
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
    assert fasta_str_to_dict(fasta) == {
        "P12345": "MALLHSARVLSGVASAFHPGLAAAASARASSWWAHVEMGPPDPILGVTEAYK",
        "P99999": "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE",
    }


def test_get_protein_sequence():
    assert (
        get_protein_sequence("P99999")
        == "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE"
    )


if "snakemake" in globals():
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        sample_peptides(
            accessions_file=snakemake.input["idmap"], output=snakemake.output[0]
        )
