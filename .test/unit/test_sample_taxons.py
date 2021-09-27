import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

# from workflow.scripts.sample_taxons import read_taxids

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_sample_taxons(tmpdir):

    # with TemporaryDirectory() as tmpdir:
    if True:
        print(tmpdir)
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".test/unit/sample_taxons/data")
        expected_path = PurePosixPath(".test/unit/sample_taxons/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("results/sample_taxons/sample_taxons.txt", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "results/sample_taxons/sample_taxons.txt",
            "-F", 
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
