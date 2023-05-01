rule convert_simulated_peptides_to_fasta:
    input:
        "results/simulation/peptides_{N}.txt",
    output:
        "results/fastas/simulated_peptides_{N}.fasta",
    log:
        "logs/simulation/convert_peptides_txt_to_fasta_{N}.txt",
    shell:
        "cat {input} | sed 's/.*/>&\\n&/' > {output} 2>{log}"
