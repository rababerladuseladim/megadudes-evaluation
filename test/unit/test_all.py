import shutil
from pathlib import Path
from test.unit import common
from test.unit.common import run_command


def test_all_dry_run(tmpdir, workflow_path):
    workdir = Path(tmpdir) / "workdir"
    data_path = Path(__file__).parent / "all"
    input_path = (data_path / "data").as_posix()
    expected_path = (data_path / "expected").as_posix()

    # Copy data to the temporary workdir.
    shutil.copytree(input_path, workdir)

    # Run the test job.
    run_command([
            "python",
            "-m",
            "snakemake",
            "-s",
            workflow_path / "workflow/Snakefile",
            "all",
            "-j1",
            "--keep-target-files",
            "--dryrun",
            "--directory",
            workdir,
        ])

    common.OutputChecker(input_path, expected_path, workdir).check()
