import random
from itertools import groupby
from pathlib import Path
from typing import cast, Dict, Iterable

from pyteomics.parser import cleave
from urllib import request


import json
import sys

LOG_HANDLE = sys.stderr


def sample_peptides(accessions_file, output):
    with open(accessions_file, "r") as handle:
        tax2acc = cast(Dict[str, list[str]], json.load(handle))

    # sample accessions
    accessions_sample = []
    for i, (tax_id, accessions) in enumerate(tax2acc.items()):
        size = min(len(accessions), 100)
        accessions_sample.extend(random.sample(accessions, size))
        if i >= 100:
            break

    print(f"Sampled {len(accessions_sample)} accessions", file=LOG_HANDLE)

    # get protein sequences for sampled accessions
    fastas = []
    chunk_size = 100
    for chunk_start in range(0, len(accessions_sample), chunk_size):
        chunk = accessions_sample[chunk_start:chunk_start + chunk_size]
        fastas.append(get_fasta(chunk))
    acc2seq = {acc: seq for f in fastas for acc, seq in convert_fasta_str_to_dict(f).items()}

    # cleave proteins sequences and sample peptides
    peptides_sample = list()
    for seq in acc2seq.values():
        peptides = cleave_protein_sequence(seq)
        peptides_sample.extend(random.sample(list(peptides), min(3, len(peptides))))

    print(f"Sampled {len(peptides_sample)} peptides", file=LOG_HANDLE)

    # write fasta
    with open(output, "w") as handle:
        handle.write("\n".join(peptides_sample))


def cleave_protein_sequence(sequence, min_length=5):
    return cleave(sequence, rule="trypsin", min_length=min_length)


def get_protein_sequence(accession):
    fasta = get_fasta([accession])
    return "".join(fasta.split("\n")[1:])


def get_fasta(uniprot_list: Iterable[str]):
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
    base_url = "https://rest.uniprot.org/"
    with request.urlopen("https://www.uniprot.org/") as f:
        uniprot_version = f.getheader('X-UniProt-Release')
    print(f"Uniprot Release Number: {uniprot_version}", file=LOG_HANDLE)
    # This makes it so we match only the ENTRY field
    uniprot_list = ["accession%3A" + id for id in uniprot_list]
    line = "+OR+".join(uniprot_list)
    url = base_url + f"uniprotkb/search?query={line}&format=fasta"
    with request.urlopen(url) as f:
        fasta = f.read().decode("utf-8").strip()
    return fasta


def convert_fasta_str_to_dict(fasta):
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
    print(tmpdir)
    sample_peptides(
        Path(__file__).parent.parent.parent / ".test/unit/sample_peptides/data/results/map_taxids_to_uniprot_accessions/tax2accessions.json",
        tmpdir / "peptides.txt"
    )


if "snakemake" in globals():
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        sample_peptides(
            accessions_file=snakemake.input["accessions"], output=snakemake.output[0]
        )
