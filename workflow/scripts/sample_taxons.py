import random
import sys
from pathlib import Path


def sample_taxons(uniprot_speclist_file, output_file, log=sys.stderr, count=10, seed=None):
    with open(uniprot_speclist_file, "r") as file_handle:
        taxids = read_taxids_from_uniprot_speclist_file(file_handle)
    if count > len(taxids):
        print(
            f"Number of requested taxids ({count}) is larger than number of available taxids ({len(taxids)}).",
            file=log,
        )
        sys.exit(1)
    random.seed(seed)
    sample_taxids = random.sample(taxids, count)
    with open(output_file, "w") as output_handle:
        for taxid in sample_taxids:
            output_handle.write(f"{taxid}\n")


def read_taxids_from_uniprot_speclist_file(file_handle):
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


def sample_taxons_multiple_times(uniprot_speclist_file, output_dir, repeats=1, log=sys.stderr):
    Path(output_dir).mkdir(exist_ok=True)
    for i in range(repeats):
        output_file = Path(output_dir, f"sample_taxons-{i}.txt")
        sample_taxons(uniprot_speclist_file, output_file, log=log)


if "snakemake" in globals():
    with open(snakemake.log[0], "w") as log_handle:
        sample_taxons_multiple_times(snakemake.input[0], snakemake.output[0], repeats=snakemake.params["repeats"], log=log_handle)
