from test.common import snakemake_run, check_output


def test_sample_peptides(tmpdir, workflow_path, prepared_workdir):
    target = "results/peptides/simulated_peptides_1.txt"
    snakefile = workflow_path / "workflow/rules/simulate_sample.smk"

    # Run the test job.
    snakemake_run(
        snakefile,
        target,
        prepared_workdir.workdir,
        additional_arguments=["--use-conda",],
    )

    check_output(prepared_workdir, mode="text")
