rule sample_peptides:
    input:
        accessions="results/map_taxids_to_uniprot_accessions/tax2accessions.json",
    output:
        "results/sample_peptides/peptides.txt",
    conda:
        "../envs/sample_peptides.yaml"
    log:
        "logs/sample_peptides.txt",
    script:
        "../scripts/sample_peptides.py"
