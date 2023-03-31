import os
import sys

import subprocess as sp
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_sample_taxons(tmpdir):
    workdir = Path(tmpdir) / "workdir"
    data_path = PurePosixPath("test/unit/sample_taxons/data")
    expected_path = PurePosixPath("test/unit/sample_taxons/expected")

    # Copy data to the temporary workdir.
    shutil.copytree(data_path, workdir)

    # dbg
    print("results/sample_taxons/sample_taxons.txt", file=sys.stderr)

    # Run the test job.
    sp.check_output([
        "python",
        "-m",
        "snakemake",
        "-s",
        "workflow/rules/sample_taxons.smk",
        "results/sample_taxons/sample_taxons.txt",
        "-j1",
        "--keep-target-files",
        "--directory",
        workdir,
    ])

    # Check the output byte by byte using cmp.
    # To modify this behavior, you can inherit from common.OutputChecker in here
    # and overwrite the method `compare_files(generated_file, expected_file),
    # also see common.py.
    common.OutputChecker(data_path, expected_path, workdir).check()