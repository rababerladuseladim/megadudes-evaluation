from test.common import snakemake_run, check_output


def test_sample_taxons(workflow_path, prepared_workdir):
    target = "results/simulation/sample_taxons_lineage_1.tsv"
    snakefile = workflow_path / "workflow/rules/simulate_sample.smk"

    # Run the test job.
    snakemake_run(
        snakefile,
        target,
        prepared_workdir.workdir
    )

    check_output(prepared_workdir, mode="text")
