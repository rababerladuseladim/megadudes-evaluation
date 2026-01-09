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
        query_fasta="results/fastas/{sample}.fasta",
        dmnd_db="results/diamond/proc/query_dbs.dmnd",
    output:
        "results/diamond/{sample}.tsv",
    log:
        "logs/diamond/run_diamond-{sample}.txt",
    conda:
        "../envs/diamond.yaml"
    threads: 256
    shell:
        # top 0: keep only hits with the same score as the top scoring one
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
        --top 0 \
        --outfmt 6 qseqid sseqid slen sstart evalue \
        -o {output} \
        > {log} 2>&1"
