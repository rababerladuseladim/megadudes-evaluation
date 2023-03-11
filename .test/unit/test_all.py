import os
import sys

import subprocess as sp
import shutil
from pathlib import PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_all_dry_run(tmpdir):
    workdir = tmpdir / "workdir"
    data_path = PurePosixPath(".test/unit/all/data")
    expected_path = PurePosixPath(".test/unit/all/expected")

    # Copy data to the temporary workdir.
    shutil.copytree(data_path, workdir)

    # Run the test job.
    sp.check_output([
        "python",
        "-m",
        "snakemake",
        "all",
        "-j1",
        "--keep-target-files",
        "--dryrun",
        "--directory",
        workdir,
    ])

    # Check the output byte by byte using cmp.
    # To modify this behavior, you can inherit from common.OutputChecker in here
    # and overwrite the method `compare_files(generated_file, expected_file),
    # also see common.py.
    common.OutputChecker(data_path, expected_path, workdir).check()
