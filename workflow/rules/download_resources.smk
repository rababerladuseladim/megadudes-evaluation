import os

rule download_uniprot_resources:
    output:
        speclist="resources/uniprot/peclist.txt",
        idmap="resources/uniprot/idmapping_selected.tab.gz",
        swissprot_fasta="resources/uniprot/swissprot.fasta.gz",
        trembl_fasta="resources/uniprot/trembl.fasta.gz",
    log:
        "logs/download_uniprot_resources.txt"
    shell:
        """\
        UNIPROT_VERSION=$(curl --head -sS https://www.uniprot.org 2>"{log}"| grep x-uniprot-release | sed 's/x-uniprot-release: //')
        echo "Uniprot Release Number: ${{UNIPROT_VERSION}}" >>"{log}"
        curl -sS https://www.uniprot.org/docs/speclist.txt > "{output.speclist}" 2>>"{log}"
        curl -sS https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz \
        > "{output.idmap}" 2>>"{log}"
        curl -sS https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz \
        > "{output.swissprot_fasta}" 2>>"{log}"
        curl -sS https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz \
        > "{output.trembl_fasta}" 2>>"{log}"
        """


rule download_ncbi_resources:
    output:
        nodes = "resources/ncbi/nodes.dmp",
        names = "resources/ncbi/names.dmp"
    log:
        "logs/download_ncbi_resources.txt"
    params:
        dir=lambda wildcards, output: os.path.dirname(output.nodes)
    shell:
        """
        curl -sS https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz 2> "{log}" | \\
        tar -xzf - --directory "{params.dir}" nodes.dmp names.dmp
        """
