import random
import sys

LOG_HANDLE = sys.stderr

def sample_lines(input_file, output_file, count=50):
    """Take a random sample of lines from the provided input file.

    The first line is expected to be the headline and written as first line to the output file.

    Args:
        input_file: path to the input file, must habe at least `count` + 1 number of lines
        output_file: path to file for writing headline and sampled lines to. Is also used for initializing
            random.seed.
        count: number of lines that should be sampled
    """
    with open(input_file, "r") as file_handle:
        headline, *lines = file_handle.readlines()
    if count > len(lines):
        print(
            f"Number of requested taxids ({count}) is larger than number of available taxids ({len(lines)}).",
            file=LOG_HANDLE,
        )
        sys.exit(1)
    random.seed(output_file)
    sampled_lines = random.sample(lines, count)
    with open(output_file, "w") as output_handle:
        output_handle.writelines([headline, *sampled_lines])


if snakemake := globals().get("snakemake"):
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        sample_lines(snakemake.input[0], snakemake.output[0])
