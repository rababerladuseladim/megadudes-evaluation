from test.common import snakemake_run, check_output


def test_map_diamond_accessions_to_taxids(tmpdir, workflow_path, prepared_workdir):
    targets = ["results/diamond/foo-bar-taxids.txt"]
    snakefile = workflow_path / "workflow/rules/qc_alignment.smk"

    # Run the test job.
    snakemake_run(
        snakefile,
        targets,
        prepared_workdir.workdir,
    )

    check_output(prepared_workdir, mode="text")
