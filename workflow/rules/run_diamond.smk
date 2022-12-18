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
        "resources/uniprot/swissprot.fasta.gz",
    output:
        "results/diamond/proc/swissprot.dmnd",
    conda:
        "../envs/diamond.yaml"
    log:
        "logs/diamond/make_diamond_db-swissprot.txt",
    threads: 256
    shell:
        "zcat {input} | diamond makedb --threads {threads} --db {output} > {log} 2>&1"


rule run_diamond:
    input:
        query_fasta="results/diamond/proc/peptides.fasta",
        dmnd_db="results/diamond/proc/swissprot.dmnd",
    output:
        "results/diamond/out.tsv",
    conda:
        "../envs/diamond.yaml"
    log:
        "logs/diamond/run_diamond.txt",
    threads: 256
    shell:
        # output fields: (for a list of available options see https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#output-options)
        # qseqid: Query Seq - id
        # sseqid: Subject Seq - id
        # slen: Subject sequence length
        # sstart: Start of alignment in subject
        # cigar: CIGAR string
        # pident: Percentage of identical matches
        # mismatch: Number of mismatches
        # evalue: Expect value
        "diamond blastp  \
        -q {input.query_fasta} \
        -d {input.dmnd_db} \
        --threads {threads}\
        --fast \
        --outfmt 6 qseqid sseqid slen sstart cigar pident mismatch evalue \
        -o {output} \
        > {log} 2>&1"
