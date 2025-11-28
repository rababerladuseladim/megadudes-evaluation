from test.common import snakemake_run, check_output


def test_mmseqs2_top_10(tmpdir, workflow_path, prepared_workdir):
    targets = ["results/mmseqs2_top_10/sample.tsv"]
    snakefile = workflow_path / "workflow/rules/run_mmseqs2.smk"

    # Run the test job.
    snakemake_run(
        snakefile,
        targets,
        prepared_workdir.workdir,
        additional_arguments=[
        "--config",
        "query_dbs=foo",]
    )

    check_output(prepared_workdir, mode="text")
