import shutil
from pathlib import Path

import pytest

from test.common import StandardPaths


@pytest.fixture
def workflow_path() -> Path:
    return Path(__file__).parent.parent


@pytest.fixture
def standard_paths(tmp_path: Path, request: pytest.FixtureRequest) -> StandardPaths:
    """Infer standard test directories based on the test file.

    The fixture expects the following directory structure next to the test file:

        data/
        expected/

    Args:
        tmp_path: A pytest-provided temporary directory path.
        request: Request context for the test.

    Returns:
        StandardPaths: Object containing workdir, input_path, and expected_path.
    """
    test_path: Path = Path(request.path)
    test_dir: Path = test_path.parent

    input_path: Path = test_dir / "data"
    expected_path: Path = test_dir / "expected"
    workdir: Path = tmp_path / "workdir"

    return StandardPaths(
        workdir=workdir,
        input_path=input_path,
        expected_path=expected_path,
    )


@pytest.fixture
def prepared_workdir(standard_paths: StandardPaths) -> StandardPaths:
    """Prepare the working directory by copying input test data into it.

    This fixture uses the ``StandardPaths`` object produced by ``standard_paths``,
    copies ``input_path`` → ``workdir``, and returns the same object so downstream
    fixtures/tests can access all paths consistently.

    Args:
        standard_paths: Standardized directory paths for the test.

    Returns:
        StandardPaths: Same object, with ``workdir`` populated with input files.
    """
    shutil.copytree(standard_paths.input_path, standard_paths.workdir)
    return standard_paths
