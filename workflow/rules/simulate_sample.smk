wildcard_constraints:
    repeat=r"\d+",


rule convert_scientific_names_to_taxonomy_ids:
    input:
        "resources/human_microbiome_species.txt",
    output:
        "results/simulation/human_microbiome_project_taxonomy_ids.txt",
    log:
        "logs/simulation/convert_scientific_names_to_taxonomy_ids.txt",
    script:
        "../scripts/convert_scientific_names_to_taxonomy_ids.py"


rule map_taxids_to_uniprot_accessions:
    input:
        idmap="resources/uniprot/idmapping_selected.tab.gz",
        tax_ids="results/simulation/human_microbiome_project_taxonomy_ids.txt",
    output:
        tax2acc_map="results/simulation/tax2accessions.json",
        tax_ids="results/simulation/human_microbiome_project_taxonomy_ids_in_uniprot.txt",
    log:
        "logs/simulation/map_taxids_to_uniprot_accessions.txt",
    resources:
        mem_gb=21,
    script:
        "../scripts/map_taxids_to_uniprot_accessions.py"


rule filter_tax_ids_and_build_lineage:
    input:
        tax_ids="results/simulation/human_microbiome_project_taxonomy_ids_in_uniprot.txt",
        ncbi_nodes="resources/ncbi/nodes.dmp",
    output:
        "results/simulation/human_microbiome_project_lineage.tsv",
    log:
        "logs/simulation/filter_tax_ids_and_build_lineage.txt",
    script:
        "../scripts/filter_tax_ids_and_build_lineage.py"


rule sample_taxons:
    input:
        "results/simulation/human_microbiome_project_lineage.tsv",
    output:
        "results/simulation/sample_taxons_lineage_{repeat}.tsv",
    log:
        "logs/simulation/sample_taxons_{repeat}.txt",
    script:
        "../scripts/sample_taxons.py"


rule sample_peptides:
    input:
        tax2acc_map="results/simulation/tax2accessions.json",
        lineage="results/simulation/sample_taxons_lineage_{repeat}.tsv",
    output:
        "results/simulation/peptides_{repeat}.txt",
    log:
        "logs/simulation/sample_peptides_{repeat}.txt",
    conda:
        "../envs/pyteomics.yaml"
    script:
        "../scripts/sample_peptides.py"


rule convert_simulated_peptides_to_fasta:
    input:
        "results/simulation/peptides_{repeat}.txt",
    output:
        "results/fastas/simulated_peptides_{repeat}.fasta",
    log:
        "logs/simulation/convert_peptides_txt_to_fasta_{repeat}.txt",
    shell:
        "cat {input} | sed 's/.*/>&\\n&/' > {output} 2>{log}"
