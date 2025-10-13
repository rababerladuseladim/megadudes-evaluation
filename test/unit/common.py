"""
Common code for unit testing of rules generated with Snakemake 6.8.0.
"""

from pathlib import Path
from typing import Literal

import os
import subprocess


def run_command(command: list[str]) -> None:
    """Run a subprocess command and print stdout/stderr if it fails."""
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        print("--- STDOUT ---")
        print(e.stdout.strip())
        print("--- STDERR ---")
        print(e.stderr.strip())
        raise
    else:
        print(result.stdout)


class OutputChecker:
    def __init__(self, data_path, expected_path, workdir, mode=Literal["binary", "text"]):
        self.data_path = data_path
        self.expected_path = expected_path
        self.workdir = workdir
        self.mode = mode

    def check(self, compare_content=True):
        """Check that all expected files are created.

        Args:
            compare_content: Compare content of expected files with created files, defaults to True
        """
        input_files = set(
            (Path(path) / f).relative_to(self.data_path)
            for path, subdirs, files in os.walk(self.data_path)
            for f in files
        )
        expected_files = set(
            (Path(path) / f).relative_to(self.expected_path)
            for path, subdirs, files in os.walk(self.expected_path)
            for f in files
        )
        unexpected_files = set()
        created_files = set()
        for path, subdirs, files in os.walk(self.workdir):
            for f in files:
                f = (Path(path) / f).relative_to(self.workdir)
                created_files.add(f)
                if str(f).startswith((".snakemake", "log")):
                    continue
                if f in expected_files:
                    if compare_content:
                        self.compare_files(self.workdir / f, self.expected_path / f)
                elif f in input_files:
                    # ignore input files
                    pass
                else:
                    unexpected_files.add(f)
        if missing_files := expected_files - created_files:
            raise ValueError(
                "Missing files:\n{}".format(
                    "\n".join(sorted(map(str, missing_files)))
                )
            )
        if unexpected_files:
            raise ValueError(
                "Unexpected files:\n{}".format(
                    "\n".join(sorted(map(str, unexpected_files)))
                )
            )

    def compare_files(self, generated_file, expected_file):

        cmd = "cmp" if self.mode == "binary" else "diff"

        try:
            subprocess.check_output([cmd, generated_file, expected_file])
        except subprocess.CalledProcessError as error:
            raise AssertionError(
                "Error comparing files:\n"
                f"    {generated_file}\n"
                f"    {expected_file}\n"
                f"Command: {' '.join(map(str, error.cmd))}\n"
                "Output:\n" +
                error.output.decode()
            )
