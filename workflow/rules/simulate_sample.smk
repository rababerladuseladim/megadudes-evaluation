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
    log:
        "logs/simulation/map_taxids_to_uniprot_accessions.txt",
    resources:
        mem_gb=21,
    script:
        "../scripts/map_taxids_to_uniprot_accessions.py"


rule convert_tax2acc_map_to_tax_ids:
    input:
        tax2acc_map="results/simulation/tax2accessions.json",
    output:
        tax_ids="results/simulation/human_microbiome_project_taxonomy_ids_in_uniprot.txt",
    run:
        import json

        with open(input.tax2acc_map) as f:
            tax2acc_map = json.load(f)
            tax_ids = list(tax2acc_map.keys())
        with open(output.tax_ids, "w") as outfile:
            outfile.write("\n".join(tax_ids))
            outfile.write("\n")


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
        "../scripts/sample_lines.py"


rule sample_peptides:
    input:
        tax2acc_map="results/simulation/tax2accessions.json",
        lineage="results/simulation/sample_taxons_lineage_{repeat}.tsv",
    output:
        "results/peptides/simulated_peptides_{repeat}_with_{percentage}_percent_noise.txt",
    log:
        "logs/simulation/sample_peptides-simulated_peptides_{repeat}_with_{percentage}_percent_noise.txt",
    conda:
        "../envs/pyteomics.yaml"
    script:
        "../scripts/sample_peptides.py"
