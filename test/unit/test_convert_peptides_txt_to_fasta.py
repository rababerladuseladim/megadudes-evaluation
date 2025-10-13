import sys
import shutil
from pathlib import Path
from test.unit import common
from test.unit.common import run_command


def test_convert_peptides_txt_to_fasta(tmpdir, workflow_path):
    workdir = Path(tmpdir) / "workdir"
    data_path = Path(__file__).parent / __name__.split(".")[-1].removeprefix("test_")
    input_path = (data_path / "data").as_posix()
    expected_path = (data_path / "expected").as_posix()

    # Copy data to the temporary workdir.
    shutil.copytree(input_path, workdir)

    # dbg
    print("results/fastas/simulated_peptides_10.fasta", file=sys.stderr)

    # Run the test job.
    command = [
        "python",
        "-m",
        "snakemake",
        "-s",
        workflow_path / "workflow/rules/convert_peptides_txt_to_fasta.smk",
        "results/fastas/simulated_peptides_10.fasta",
        "-j1",
        "--keep-target-files",
        "--directory",
        workdir,
    ]
    run_command(command)
    common.OutputChecker(input_path, expected_path, workdir).check()
