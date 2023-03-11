import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_map_taxids_to_uniprot_accessions():
    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".test/unit/map_taxids_to_uniprot_accessions/data")
        expected_path = PurePosixPath(".test/unit/map_taxids_to_uniprot_accessions/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("results/map_taxids_to_uniprot_accessions/tax2accessions.json", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake",
            "-s",
            "workflow/rules/map_taxids_to_uniprot_accessions.smk",
            "results/map_taxids_to_uniprot_accessions/tax2accessions.json",
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
