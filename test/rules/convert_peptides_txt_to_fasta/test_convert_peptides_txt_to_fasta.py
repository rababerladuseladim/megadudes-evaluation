from test.common import snakemake_run, check_output


def test_convert_peptides_txt_to_fasta(tmpdir, workflow_path, prepared_workdir):
    target = "results/fastas/simulated_peptides_10.fasta"
    snakefile = workflow_path / "workflow/rules/convert_peptides_txt_to_fasta.smk"

    # Run the test job.
    snakemake_run(
        snakefile,
        target,
        prepared_workdir.workdir,
    )

    check_output(prepared_workdir)
