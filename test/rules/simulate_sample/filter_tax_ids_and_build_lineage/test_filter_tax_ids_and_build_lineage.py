from test.common import snakemake_run, check_output

def test_filter_tax_ids_and_build_lineage(tmpdir, workflow_path, prepared_workdir):
    target = "results/simulation/human_microbiome_project_lineage.tsv"
    snakefile = workflow_path / "workflow/rules/simulate_sample.smk"

    # Run the test job.
    snakemake_run(
        snakefile,
        target,
        prepared_workdir.workdir,
    )

    check_output(prepared_workdir, mode="text")
