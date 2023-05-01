rule convert_msfragger_peptides_tsv_to_fasta:
    input:
        config["msfragger_peptides_tsv"],
    output:
        "results/diamond/proc/peptides.fasta",
    log:
        "logs/diamond/convert_msfragger_peptides_tsv_to_fasta.txt",
    shell:
        "cat {input} | tail -n +2 | cut -f1 | sed 's/.*/>&\\n&/' > {output} 2>{log}"
