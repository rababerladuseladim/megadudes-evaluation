rule convert_msfragger_peptides_tsv_to_fasta:
    input:
        lambda wc: samples.loc[wc.sample_name, "msfragger_peptides_tsv"],
    output:
        "results/fastas/sample_{sample_name}.fasta",
    log:
        "logs/fastas/convert_msfragger_peptides_tsv_to_fasta-sample_{sample_name}.txt",
    shell:
        "cat {input} | tail -n +2 | cut -f1 | sed 's/.*/>&\\n&/' > {output} 2>{log}"
