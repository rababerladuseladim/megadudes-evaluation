rule sample_taxons:
    input:
        "resources/uniprot/speclist.txt",
    output:
        "results/simulation/sample_taxons_{repeat}.txt",
    log:
        "logs/simulation/sample_taxons_{repeat}.txt",
    script:
        "../scripts/sample_taxons.py"


rule map_taxids_to_uniprot_accessions:
    input:
        idmap="resources/uniprot/idmapping_selected.tab.gz",
        taxids="results/simulation/sample_taxons_{repeat}.txt",
    output:
        "results/simulation/tax2accessions_{repeat}.json",
    log:
        "logs/simulation/map_taxids_to_uniprot_accessions_{repeat}.txt",
    script:
        "../scripts/map_taxids_to_uniprot_accessions.py"


rule sample_peptides:
    input:
        accessions="results/simulation/tax2accessions_{repeat}.json",
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