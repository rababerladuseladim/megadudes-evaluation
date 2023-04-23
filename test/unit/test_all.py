import subprocess as sp
import shutil
from pathlib import Path
from test.unit import common


def test_all_dry_run(tmpdir):
    workdir = Path(tmpdir) / "workdir"
    data_path = Path(__file__).parent / "all"
    input_path = (data_path / "data").as_posix()
    expected_path = (data_path / "expected").as_posix()

    # Copy data to the temporary workdir.
    shutil.copytree(input_path, workdir)

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
    common.OutputChecker(input_path, expected_path, workdir).check()
