rule convert_msfragger_peptides_tsv_to_fasta:
    input:
        config["msfragger_peptides_tsv"],
    output:
        "results/diamond/proc/peptides.fasta",
    log:
        "logs/diamond/convert_msfragger_peptides_tsv_to_fasta.txt",
    shell:
        "cat {input} | tail -n +2 | cut -f1 | sed 's/.*/>&\\n&/' > {output} 2>{log}"


rule make_diamond_db:
    input:
        config["query_dbs"],
    output:
        "results/diamond/proc/query_dbs.dmnd",
    log:
        "logs/diamond/make_diamond_db-query_dbs.txt",
    conda:
        "../envs/diamond.yaml"
    threads: 256
    shell:
        "zcat {input} | diamond makedb --threads {threads} --db {output} > {log} 2>&1"


rule run_diamond:
    input:
        query_fasta="results/diamond/proc/peptides.fasta",
        dmnd_db="results/diamond/proc/query_dbs.dmnd",
    output:
        "results/diamond/out.tsv",
    log:
        "logs/diamond/run_diamond.txt",
    conda:
        "../envs/diamond.yaml"
    threads: 256
    shell:
        # output fields: (for a list of available options see https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#output-options)
        # qseqid: Query Seq - id
        # sseqid: Subject Seq - id
        # slen: Subject sequence length
        # sstart: Start of alignment in subject
        # evalue: Expect value
        "diamond blastp  \
        -q {input.query_fasta} \
        -d {input.dmnd_db} \
        --threads {threads}\
        --fast \
        --outfmt 6 qseqid sseqid slen sstart evalue \
        -o {output} \
        > {log} 2>&1"
