from test.common import snakemake_run, check_output


def test_build_diamond_tax_id_lineage(tmpdir, workflow_path, prepared_workdir):
    targets = ["results/diamond/foo-bar-lineage.tsv"]
    snakefile = workflow_path / "workflow/rules/qc_alignment.smk"
    configfile = prepared_workdir.workdir / "config" / "config.yaml"

    # Run the test job.
    snakemake_run(
        snakefile,
        targets,
        prepared_workdir.workdir,
        additional_arguments=[
            "--configfile",
            configfile.as_posix(),
        ]
    )

    check_output(prepared_workdir, mode="text")
