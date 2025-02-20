import sys

import subprocess as sp
import shutil
from pathlib import Path
from test.unit import common


def test_sample_peptides(tmpdir, workflow_path):
    target = "results/simulation/peptides_1.txt"

    workdir = Path(tmpdir) / "workdir"
    print("workdir:", workdir)
    data_path = Path(__file__).parent / __name__.split(".")[-1].removeprefix("test_")
    input_path = (data_path / "data").as_posix()
    expected_path = (data_path / "expected").as_posix()

    # Copy data to the temporary workdir.
    shutil.copytree(input_path, workdir)

    # dbg
    print(target, file=sys.stderr)

    # Run the test job.
    sp.check_output([
        "python",
        "-m",
        "snakemake",
        "--use-conda",
        "-s",
        workflow_path / "workflow/rules/simulate_sample.smk",
        target,
        "-j1",
        "--directory",
        workdir,
    ])

    # Check the output byte by byte using cmp.
    # To modify this behavior, you can inherit from common.OutputChecker in here
    # and overwrite the method `compare_files(generated_file, expected_file),
    # also see common.py.
    common.OutputChecker(input_path, expected_path, workdir, mode="text").check(compare_content=True)
