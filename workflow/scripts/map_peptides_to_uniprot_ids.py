from collections import defaultdict
import json
import urllib.parse
import urllib.request
import sys

LOG_HANDLE = sys.stderr


def map_peptides_to_uniprot_ids(peptide_file, out_path):
    # get uniprot database version
    with urllib.request.urlopen('https://unipept.ugent.be/private_api/metadata') as r:
        data = r.read()
    unipept_resp = json.loads(data)
    db_version = unipept_resp["db_version"]
    print(f"unipept database version: {db_version}", file=LOG_HANDLE)

    # map peptides
    url = 'http://api.unipept.ugent.be/api/v1/pept2prot.json'
    trypticPeptides = []
    # todo: limit to 100 peptide per api call
    # from unipept doku: When performing bulk searches, we suggest splitting the input set over requests of 100 peptides each.
    with open(peptide_file) as infile:
        for line in infile.readlines():
            content = line.strip()
            if content:
                trypticPeptides.append(content)

    post_data = [('equate_il', 'true')]
    post_data += [('input[]', x) for x in trypticPeptides]
    data = urllib.parse.urlencode(post_data).encode()
    req = urllib.request.Request(url, data=data, method='POST')

    with urllib.request.urlopen(req) as response:
        rdata = response.read()

    unipept_resp = json.loads(rdata)

    pep2prot = defaultdict(list)
    for res in unipept_resp:
        pep2prot[res["peptide"]].append(res["uniprot_id"])

    with open(out_path, "w") as outfile:
        json.dump(pep2prot, outfile)


if "snakemake" in globals():
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        map_peptides_to_uniprot_ids(snakemake.input[0], snakemake.output[0])
