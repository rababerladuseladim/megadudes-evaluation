rule download_uniprot_idmappig:
    output:
        "uniprot/idmapping_selected.tab.gz"
    log:
        "logs/uniprot/idmapping_download.txt"
    shell:
        """\
        wget \
        -o \"{log}\" \
        -O \"{output[0]}\" \
        https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz"""

rule map_taxids_to_uniprot_accessions:
    input:
        idmap="uniprot/idmapping_selected.tab.gz",
        taxids="results/sample_taxons/sample_taxons.txt"
    output:
        "results/map_taxids_to_uniprot_accessions/tax2accessions.json"
    log:
        "logs/map_taxids_to_uniprot_accssions.txt"
    script:
        "../scripts/map_taxids_to_uniprot_accessions.py"
