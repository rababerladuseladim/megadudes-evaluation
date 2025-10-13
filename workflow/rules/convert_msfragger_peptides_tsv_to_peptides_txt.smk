rule convert_msfragger_peptides_tsv_to_peptides_txt:
    input:
        lambda wc: samples.loc[wc.sample_name, "msfragger_peptides_tsv"],
    output:
        "results/peptides/sample_{sample_name}.txt",
    log:
        "logs/peptides/convert_msfragger_peptides_tsv_to_plain_peptides_txt-sample_{sample_name}.txt",
    shell:
        "cat {input} | tail -n +2 | cut -f1 > {output} 2>{log}"
