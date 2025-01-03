import random
import sys


def sample_taxons(input_file, output_file, log=sys.stderr, count=100):
    """Take a random sample of lines from the provided input file.

    The first line is expected to be the headline and written as first line to the output file.

    Args:
        input_file: path to the input file, must habe at least `count` + 1 number of lines
        output_file: path to file for writing TaxIDs to, one TaxID per line. Is also used for initializing
            random.seed.
        log: handle for writing the log
        count: number of taxa that should be sampled
    """
    with open(input_file, "r") as file_handle:
        headline, *lines = file_handle.readlines()
    if count > len(lines):
        print(
            f"Number of requested taxids ({count}) is larger than number of available taxids ({len(lines)}).",
            file=log,
        )
        sys.exit(1)
    random.seed(output_file)
    sample_lines = random.sample(lines, count)
    with open(output_file, "w") as output_handle:
        output_handle.writelines([headline, *sample_lines])


if snakemake := globals().get("snakemake"):
    with open(snakemake.log[0], "w") as log_handle:
        sample_taxons(snakemake.input[0], snakemake.output[0], log=log_handle)
