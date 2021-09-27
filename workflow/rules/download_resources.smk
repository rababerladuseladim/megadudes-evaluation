rule download_uniprot_idmapping:
    output:
        idmap="resources/idmapping_selected.tab.gz"
    log:
        "logs/download_uniprot_idmapping.txt"
    shell:
        """\
        wget \
        -o \"{log}\" \
        -O \"{output.idmap}\" \
        https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz"""

rule download_uniprot_speclist:
    output:
        speclist="resources/speclist.txt"
    log:
        "logs/download_uniprot_speclist.txt"
    shell:
        """\
        wget \
        -o \"{log}\" \
        -O \"{output.speclist}\" \
        https://www.uniprot.org/docs/speclist.txt"""

