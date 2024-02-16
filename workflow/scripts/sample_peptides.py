import random
import numpy as np
from itertools import groupby
from typing import cast, Dict, Iterable

from numpy.random import Generator as NumpyGenerator
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
        chunk = accessions_sample[chunk_start : chunk_start + chunk_size]
        fastas.append(connector.get_fasta(chunk))
    acc2seq = {
        acc: seq for f in fastas for acc, seq in convert_fasta_str_to_dict(f).items()
    }

    # ensure deterministic behaviour although previous fetching of sequences takes a non-deterministic amount of
    # computation
    random.setstate(seed)
    numpy_random_number_generator = np.random.default_rng(seed=random.getstate()[1][0])

    # cleave proteins sequences and sample peptides
    sequences = sorted(acc2seq.values())
    peptides_sample = list()
    for seq in sequences:
        peptides_sample.extend(
            sample_peptides_from_sequence(numpy_random_number_generator, seq)
        )

    print(f"Sampled {len(peptides_sample)} peptides", file=LOG_HANDLE)

    # write fasta
    with open(output, "w") as handle:
        handle.write("\n".join(peptides_sample) + "\n")


def sample_peptides_from_sequence(
    numpy_random_number_generator: NumpyGenerator, sequence: str
) -> list[str]:
    """Digest sequence according to tryptic digestion rules and return a random sample of unique peptides.

    The sample size is drawn from a poisson distribution with mean 10.
    Digestion allows up to 2 missed cleavages, with the weights of each missed cleavage as follows:
    - 0 missed cleavages: 100
    - 1: 10
    - 2: 1

    Args:
        numpy_random_number_generator: initialized numpy Generator
        sequence: amino acid sequence to be digested

    Returns:
        sample of peptides
    """
    (
        peptides_0_missed_cleavages,
        peptides_1_missed_cleavages,
        peptides_2_missed_cleavages,
    ) = [
        list(
            sorted(cleave_protein_sequence(sequence, missed_cleavages=missed_cleavages))
        )
        for missed_cleavages in [0, 1, 2]
    ]
    sample_size = numpy_random_number_generator.poisson(lam=10)
    missed_cleavages_distribution = random.choices(
        [0, 1, 2], weights=[100, 10, 1], k=sample_size
    )

    return (
        random.sample(
            population=peptides_0_missed_cleavages,
            k=min(
                missed_cleavages_distribution.count(0), len(peptides_0_missed_cleavages)
            ),
        )
        + random.sample(
            population=peptides_1_missed_cleavages,
            k=min(
                missed_cleavages_distribution.count(1), len(peptides_1_missed_cleavages)
            ),
        )
        + random.sample(
            population=peptides_2_missed_cleavages,
            k=min(
                missed_cleavages_distribution.count(2), len(peptides_2_missed_cleavages)
            ),
        )
    )


def cleave_protein_sequence(
    sequence: str,
    min_length: int = 7,
    max_length: int = 30,
    missed_cleavages: int = 0,
    *args,
    **kwargs,
) -> set[str]:
    kwargs.update(
        {
            "sequence": sequence,
            "rule": "trypsin",
            "min_length": min_length,
            "max_length": max_length,
        }
    )
    if missed_cleavages > 0:
        return cleave(*args, missed_cleavages=missed_cleavages, **kwargs) - cleave(
            *args, missed_cleavages=missed_cleavages - 1, **kwargs
        )
    return cleave(*args, missed_cleavages=missed_cleavages, **kwargs)


def get_protein_sequence(accession) -> str:
    connector = UniProtConnector()
    fasta = connector.get_fasta([accession])
    return "".join(fasta.split("\n")[1:])


class UniProtConnector:
    url = "https://rest.uniprot.org/uniprotkb/"

    def __init__(self):
        self.session = requests.Session()
        self.uniprot_version = requests.get(self.url).headers.get("X-UniProt-Release")
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


if snakemake := globals().get("snakemake"):
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        sample_peptides(
            accessions_file=snakemake.input["accessions"], output=snakemake.output[0]
        )
