import shutil
import sys

from pathlib import Path
from test.unit import common
from test.unit.common import run_command


def test_sample_peptides(tmpdir, workflow_path):
    target = "results/peptides/simulated_peptides_1.txt"

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
    run_command([
        "python",
        "-m",
        "snakemake",
        "--use-conda",
        "-s",
        workflow_path / "workflow/rules/simulate_sample.smk",
        target,
        "-j1",
        "--keep-target-files",
        "--directory",
        workdir,
    ])

    common.OutputChecker(input_path, expected_path, workdir, mode="text").check(compare_content=True)
