from test.common import snakemake_run, check_output


def test_all_dry_run(tmpdir, workflow_path, prepared_workdir):

    target = "all"
    snakefile = workflow_path / "workflow/Snakefile"

    # Run the test job.
    snakemake_run(
        snakefile,
        target,
        prepared_workdir.workdir,
        additional_arguments=["--dryrun"]
    )

    check_output(prepared_workdir)
