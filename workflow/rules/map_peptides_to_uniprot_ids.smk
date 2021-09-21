# An example collection of Snakemake rules imported in the main Snakefile.
rule map_peptides_to_uniprot_ids:
    input:
        config["peptides"]
    output:
        "results/map_peptides_to_uniprot_ids/pep2uniprot_ids.json"
    log:
        "results/map_peptides_to_uniprot_ids/log.txt"
    script:
        "../scripts/map_peptides_to_uniprot_ids.py"
