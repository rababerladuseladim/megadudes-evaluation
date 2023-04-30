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
