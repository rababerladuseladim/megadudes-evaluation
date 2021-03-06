rule convert_msfragger_peptides_tsv_to_fasta:
    input: config["msfragger_peptides_tsv"]
    output: "results/diamond/proc/peptides.fasta"
    log: "logs/diamond/convert_msfragger_peptides_tsv_to_fasta.txt"
    shell: "cat {input} | tail -n +2 | cut -f1 | sed 's/.*/>&\\n&/' > {output} 2>{log}"


rule make_diamond_db:
    input: "resources/uniprot/swissprot.fasta.gz"
    output: "results/diamond/proc/swissprot.dmnd"
    log: "logs/diamond/make_diamond_db-swissprot.txt"
    conda: "../envs/diamond.yaml"
    shell: "zcat {input} | diamond makedb --db {output} > {log} 2>&1"


rule run_diamond:
    input:
        query_fasta="results/diamond/proc/peptides.fasta",
        dmnd_db="results/diamond/proc/swissprot.dmnd"
    output: "results/diamond/out.tsv"
    conda: "../envs/diamond.yaml"
    log: "logs/diamond/run_diamond.txt"
    shell:
      "diamond blastp  -q {input.query_fasta} -d {input.dmnd_db} --fast --outfmt 6 qseqid sseqid slen sstart cigar pident -o {output} > {log} 2>&1"


rule plot_histogram_of_hits:
    input: "results/diamond/out.tsv"
    output: "plots/diamond/hist_of_hits.svg"
    log: "logs/diamond/plot_histogram_of_hits.txt"
    conda:
        "../envs/qc_plots.yaml"
    script:
        "../scripts/diamond-hist_of_hits.py"

