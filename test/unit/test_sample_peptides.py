import os
import sys

import subprocess as sp
import shutil
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_sample_peptides(tmpdir):
    workdir = Path(tmpdir) / "workdir"
    print("workdir:", workdir)
    data_path = Path(__file__).parent / __name__.removeprefix("test_")
    input_path = (data_path / "data").as_posix()
    expected_path = (data_path / "expected").as_posix()

    # Copy data to the temporary workdir.
    shutil.copytree(input_path, workdir)

    # dbg
    print("results/sample_peptides/peptides.txt", file=sys.stderr)

    # Run the test job.
    sp.check_output([
        "python",
        "-m",
        "snakemake",
        "--use-conda",
        "-s",
        common.WORKFLOW_PATH / "workflow/rules/sample_peptides.smk",
        "results/sample_peptides/peptides.txt",
        "-j1",
        "--keep-target-files",
        "--directory",
        workdir,
    ])

    # Check the output byte by byte using cmp.
    # To modify this behavior, you can inherit from common.OutputChecker in here
    # and overwrite the method `compare_files(generated_file, expected_file),
    # also see common.py.
    common.OutputChecker(input_path, expected_path, workdir).check(compare_content=False)
