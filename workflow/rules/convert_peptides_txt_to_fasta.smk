rule convert_peptides_txt_to_fasta:
    input:
        "results/peptides/{sample_name}.txt",
    output:
        "results/fastas/{sample_name}.fasta",
    log:
        "logs/fastas/convert_peptides_txt_to_fasta-{sample_name}.txt",
    shell:
        "cat {input} | sed 's/.*/>&\\n&/' > {output} 2>{log}"
