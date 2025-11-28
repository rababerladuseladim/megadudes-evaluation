"""
Common code for unit testing of rules generated with Snakemake 6.8.0.
"""
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Iterable

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
    def __init__(self, data_path, expected_path, workdir, mode: Literal["binary", "text"] = "text"):
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
            ) from None


@dataclass
class StandardPaths:
    """Container for standardized test directory paths.

    Attributes:
        workdir: Temporary working directory where Snakemake runs.
        input_path: Path to the test's input data.
        expected_path: Path to the expected output data.
    """
    workdir: Path
    input_path: Path
    expected_path: Path


def check_output(standard_paths: StandardPaths, mode: Literal["binary", "text"] = "text"):
    OutputChecker(standard_paths.input_path, standard_paths.expected_path, standard_paths.workdir,
                  mode=mode).check()


def snakemake_run(
    snakefile: Path,
    target: str,
    workdir: Path,
    additional_arguments: Iterable[str] = None,
) -> None:
    """Run Snakemake with a given workflow, target, and working directory.

    Args:
        snakefile: Path to the Snakemake workflow (.smk) file.
        target: The Snakemake rule target to build.
        workdir: Directory in which Snakemake will be executed.
        additional_arguments: Optional extra command-line
            arguments to append to the Snakemake invocation. Useful for passing
            flags such as ``--config`` or ``--touch``. Defaults to ``None``.

    Example:
        snakemake_run(
            snakefile=workflow_path / "rules/a.smk",
            target="results/a.txt",
            workdir=workdir,
            additional_arguments=["--config", "x=1"]
        )
    """
    args = [
        "python",
        "-m",
        "snakemake",
        "-s",
        snakefile,
        target,
        "-j1",
        "--keep-target-files",
        "--directory",
        workdir,
    ]

    if additional_arguments:
        args.extend(additional_arguments)

    run_command(args)
