import os.path

rule download_uniprot_idmapping:
    output:
        idmap="resources/idmapping_selected.tab.gz"
    log:
        "logs/download_uniprot_idmapping.txt"
    shell:
        """\
        UNIPROT_VERSION=$(curl --head -sS https://www.uniprot.org 2>"{log}"| grep x-uniprot-release | sed 's/x-uniprot-release: //')
        echo "Uniprot Release Number: ${{UNIPROT_VERSION}}" >>"{log}"
        curl -sS https://www.uniprot.org/docs/speclist.txt > "{output.idmap}" 2>>"{log}"
        """

rule download_uniprot_speclist:
    output:
        speclist="resources/speclist.txt"
    log:
        "logs/download_uniprot_speclist.txt"
    shell:
        """\
        UNIPROT_VERSION=$(curl --head -sS https://www.uniprot.org 2>"{log}"| grep x-uniprot-release | sed 's/x-uniprot-release: //')
        echo "Uniprot Release Number: ${{UNIPROT_VERSION}}" >>"{log}"
        curl -sS https://www.uniprot.org/docs/speclist.txt > "{output.speclist}" 2>>"{log}"
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
