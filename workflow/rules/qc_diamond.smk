rule plot_diamond_qc:
    input:
        "results/diamond/{sample}.tsv",
        "results/fastas/{sample}.fasta",
    output:
        report("plots/diamond/qc-{sample}.svg"),
    log:
        "logs/diamond/plot_diamond_qc-{sample}.txt",
    conda:
        "../envs/qc_plots.yaml"
    script:
        "../scripts/plot_diamond_qc.py"


rule extract_diamond_accessions:
    input:
        "results/diamond/{sample}.tsv",
    output:
        "results/diamond/{sample}-accessions.txt",
    log:
        "logs/diamond/extract_diamond_accessions-{sample}.txt",
    shell:
        "cut -f 2 {input} | grep -o '|.*|' | grep -Eo '[^|]*' > {output} 2>{log}"


rule map_diamond_accessions_to_taxids:
    input:
        idmap="resources/uniprot/idmapping_selected.tab.gz",
        accessions="results/diamond/{sample}-accessions.txt",
    output:
        "results/diamond/{sample}-taxids.txt",
    log:
        "logs/diamond/map_uniprot_accessions_to_taxids-{sample}.txt",
    resources:
        mem_gb=21,
    script:
        "../scripts/map_uniprot_accessions_to_taxids.py"


rule build_diamond_tax_id_lineage:
    input:
        tax_ids="results/diamond/{sample}-taxids.txt",
        ncbi_nodes="resources/ncbi/nodes.dmp",
    output:
        lineage="results/diamond/{sample}-lineage.tsv",
    log:
        "logs/diamond/build_diamond_tax_id_lineage_{sample}.txt",
    script:
        "../scripts/filter_tax_ids_and_build_lineage.py"
