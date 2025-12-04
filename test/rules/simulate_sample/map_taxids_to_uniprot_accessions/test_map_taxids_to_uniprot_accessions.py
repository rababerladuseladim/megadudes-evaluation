from test.common import snakemake_run, check_output


def test_map_taxids_to_uniprot_accessions(tmpdir, workflow_path, prepared_workdir):
    targets = ["results/simulation/tax2accessions.json"]
    snakefile = workflow_path / "workflow/rules/simulate_sample.smk"

    # Run the test job.
    snakemake_run(
        snakefile,
        targets,
        prepared_workdir.workdir,
    )

    check_output(prepared_workdir, mode="text")
