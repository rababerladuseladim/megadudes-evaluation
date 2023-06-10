rule sample_taxons:
    input:
        "resources/uniprot/speclist.txt",
    output:
        "results/sample_taxons/sample_taxons.txt",
    log:
        "logs/sample_taxons.txt",
    script:
        "../scripts/sample_taxons.py"


rule map_taxids_to_uniprot_accessions:
    input:
        idmap="resources/uniprot/idmapping_selected.tab.gz",
        taxids="results/sample_taxons/sample_taxons.txt",
    output:
        "results/map_taxids_to_uniprot_accessions/tax2accessions.json",
    log:
        "logs/map_taxids_to_uniprot_accssions.txt",
    script:
        "../scripts/map_taxids_to_uniprot_accessions.py"


rule sample_peptides:
    input:
        accessions="results/map_taxids_to_uniprot_accessions/tax2accessions.json",
    output:
        "results/sample_peptides/peptides.txt",
    log:
        "logs/sample_peptides.txt",
    conda:
        "../envs/pyteomics.yaml"
    script:
        "../scripts/sample_peptides.py"


rule convert_simulated_peptides_to_fasta:
    input:
        "results/simulation/peptides_{N}.txt",
    output:
        "results/fastas/simulated_peptides_{N}.fasta",
    log:
        "logs/simulation/convert_peptides_txt_to_fasta_{N}.txt",
    shell:
        "cat {input} | sed 's/.*/>&\\n&/' > {output} 2>{log}"
