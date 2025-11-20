import sys
import shutil
from pathlib import Path
from test.unit import common
from test.unit.common import run_command


def test_convert_scientific_names_to_taxonomy_ids(tmpdir, workflow_path):
    target = "results/simulation/human_microbiome_project_taxonomy_ids.txt"

    workdir = Path(tmpdir) / "workdir"
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

    # Check the output byte by byte using cmp.
    # To modify this behavior, you can inherit from common.OutputChecker in here
    # and overwrite the method `compare_files(generated_file, expected_file),
    # also see common.py.
    common.OutputChecker(input_path, expected_path, workdir, mode="text").check(compare_content=True)
