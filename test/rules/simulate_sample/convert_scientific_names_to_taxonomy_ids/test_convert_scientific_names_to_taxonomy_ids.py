from test.common import snakemake_run, check_output

def test_convert_scientific_names_to_taxonomy_ids(tmpdir, workflow_path, prepared_workdir):
    targets = ["results/simulation/human_microbiome_project_taxonomy_ids.txt"]
    snakefile = workflow_path / "workflow/rules/simulate_sample.smk"

    # Run the test job.
    snakemake_run(
        snakefile,
        targets,
        prepared_workdir.workdir,
        additional_arguments=[
        "--use-conda",]
    )

    check_output(prepared_workdir, mode="text")
