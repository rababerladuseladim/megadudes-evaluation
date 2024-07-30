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


rule create_mmseqs2_index:
    input:
        mmseqs2_result="results/mmseqs2/{sample}.tsv",
    output:
        "results/mmseqs2/{sample}.tsv.index",
    log:
        "logs/mmseqs2/create_mmseqs2_index-{sample}.txt",
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        'mmseqs createindex {input.mmseqs2_result} "$TMPDIR" > {log} 2>&1'


rule filter_mmseqs2:
    input:
        mmseqs2_result="results/mmseqs2/{sample}.tsv",
        index="results/mmseqs2/{sample}.tsv.index",
    output:
        "results/mmseqs2/filtered/{sample}.tsv",
    log:
        "logs/mmseqs2/filter_mmseqs2-{sample}.txt",
    conda:
        "../envs/mmseqs2.yaml"
    threads: 256
    shell:
        "mmseqs filterdb \
        {input.mmseqs2_result} \
        {output} \
        --extract-lines 10 \
        > {log} 2>&1"
