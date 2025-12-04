from test.common import snakemake_run, check_output


def test_sample_peptides(tmpdir, workflow_path, prepared_workdir):
    targets = [
        "results/peptides/simulated_peptides_1_with_0_percent_noise.txt",
        "results/peptides/simulated_peptides_1_with_1_percent_noise.txt",
    ]
    snakefile = workflow_path / "workflow/rules/simulate_sample.smk"
    conda_base_path = workflow_path / ".snakemake" / "conda"

    # Run the test job.
    snakemake_run(
        snakefile,
        targets,
        prepared_workdir.workdir,
        additional_arguments=[
            "--use-conda",
            "--conda-base-path",
            conda_base_path.as_posix(),
        ],
    )

    check_output(prepared_workdir, mode="text")
