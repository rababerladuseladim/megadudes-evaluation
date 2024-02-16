rule make_mmseqs2_db:
    input:
        config["query_dbs"],
    output:
        "results/mmseqs2/proc/query_dbs",
    log:
        "logs/mmseqs2/make_mmseqs2_db-query_dbs.txt",
    conda:
        "../envs/mmseqs2.yaml"
    threads: 256
    shell:
        """
        MMSEQS_NUM_THREADS={threads};
        mmseqs createdb {input} {output} > {log} 2>&1;
        mmseqs createindex {output} "$TMPDIR"
        """


rule run_mmseqs2:
    input:
        query_fasta="results/fastas/{sample}.fasta",
        mmseqs2_db="results/mmseqs2/proc/query_dbs",
    output:
        "results/mmseqs2/{sample}.tsv",
    log:
        "logs/mmseqs2/run_mmseqs2-{sample}.txt",
    conda:
        "../envs/mmseqs2.yaml"
    threads: 256
    shell:
        'mmseqs easy-search \
        {input.query_fasta} \
        {input.mmseqs2_db} \
        {output} \
        "$TMPDIR" \
        --threads {threads} \
        --seed-sub-mat VTML40.out \
        -s 2 \
        --comp-bias-corr 0 \
        --mask 0 \
        --gap-open 16 \
        --gap-extend 2 \
        --min-length 9 \
        --spaced-kmer-pattern 11011101 \
        -e 10000000 \
        -k 6 \
        --format-output "query,target,tlen,tstart,evalue" \
        > {log} 2>&1'
