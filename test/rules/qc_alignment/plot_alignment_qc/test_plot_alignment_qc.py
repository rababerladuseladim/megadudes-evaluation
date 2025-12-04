from test.common import snakemake_run, check_output


def test_plot_alignment_qc(tmpdir, workflow_path, prepared_workdir):
    targets = ["plots/alignment/qc-foo.svg"]
    snakefile = workflow_path / "workflow/rules/qc_alignment.smk"
    configfile = prepared_workdir.workdir / "config" / "config.yaml"
    conda_base_path = workflow_path / ".snakemake" / "conda"

    # Run the test job.
    snakemake_run(
        snakefile,
        targets,
        prepared_workdir.workdir,
        additional_arguments=[
            "--use-conda",
            "--configfile",
            configfile.as_posix(),
            "--conda-base-path",
            conda_base_path.as_posix(),
        ]
    )

    check_output(prepared_workdir, mode="presence")
