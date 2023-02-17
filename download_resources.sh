#!/bin/bash

# set bash to strict mode http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

# uniprot download
mkdir -p resources/uniprot
UNIPROT_VERSION=$(curl --head -sS https://rest.uniprot.org/uniprotkb/ | grep x-uniprot-release: | sed 's/x-uniprot-release: //');
echo "Uniprot Release Number: $UNIPROT_VERSION" > resources/uniprot/version.txt;
curl --output-dir resources/uniprot \
    --remote-name https://www.uniprot.org/docs/speclist.txt \
    --remote-name https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz \
    --remote-name https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz \
    --remote-name https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz;

# ncbi download
mkdir -p resources/ncbi
echo "NCBI data downloaded on: $(date)" > resources/ncbi/version.txt;
curl https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | \
        tar -xzf - --directory resources/ncbi nodes.dmp names.dmp
