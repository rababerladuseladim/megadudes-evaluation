rule map_peptides_to_uniprot_ids:
    input:
        config["peptides"]
    output:
        "results/map_peptides_to_uniprot_ids/pep2uniprot_ids.json"
    log:
        "logs/map_peptides_to_uniprot_ids.txt"
    script:
        "../scripts/map_peptides_to_uniprot_ids.py"
