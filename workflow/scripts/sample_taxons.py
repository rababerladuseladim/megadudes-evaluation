import random
import sys
from typing import Literal


def sample_taxons(input_file, output_file, log=sys.stderr, file_type: Literal["uniprot_speclist", "plain"] = "uniprot_speclist", count=100):
    """Take a random sample of taxids from the uniprot speclist file.

    Args:
        input_file: path to file containing one taxid per line or uniprot speclist file, get from: https://www.uniprot.org/docs/speclist.txt
        output_file: path to file for writing TaxIDs to, one TaxID per line. Is also used for initializing
            random.seed.
        log: handle for writing the log
        file_type: type of the input_file
        count: number of taxa that should be sampled
    """
    with open(input_file, "r") as file_handle:
        if file_type == "uniprot_speclist":
            taxids = read_taxids_from_uniprot_speclist_file(file_handle)
        else:
            taxids = [l.strip() for l in file_handle]
    if count > len(taxids):
        print(
            f"Number of requested taxids ({count}) is larger than number of available taxids ({len(taxids)}).",
            file=log,
        )
        sys.exit(1)
    random.seed(output_file)
    sample_taxids = random.sample(taxids, count)
    with open(output_file, "w") as output_handle:
        for taxid in sample_taxids:
            output_handle.write(f"{taxid}\n")


def read_taxids_from_uniprot_speclist_file(file_handle) -> list[int]:
    """

    :param file_handle: open file handle to uniprot speclist file. download from: https://www.uniprot.org/docs/speclist.txt
    :return:
    """

    taxids = []
    is_headline_reached = False
    is_taxid_section = False
    for line in file_handle.readlines():
        if not is_taxid_section:
            if not is_headline_reached:
                if line.startswith("(1) Real organism codes"):
                    is_headline_reached = True
            else:
                if line.startswith("_____"):
                    is_taxid_section = True
            continue
        else:
            # stop at empty line
            if not line.strip():
                break
            line_up_to_taxid = line[:15]
            if line_up_to_taxid.strip():
                taxid = int(line_up_to_taxid.split(" ")[-1])
                taxids.append(taxid)

    return taxids


if snakemake := globals().get("snakemake"):
    with open(snakemake.log[0], "w") as log_handle:
        sample_taxons(snakemake.input[0], snakemake.output[0], file_type="plain", log=log_handle)
