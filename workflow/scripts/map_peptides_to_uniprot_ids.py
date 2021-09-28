from collections import defaultdict
import json
import urllib.parse
import urllib.request
import sys

LOG_HANDLE = sys.stderr


def map_peptides_to_uniprot_ids(peptide_file, out_path):
    # get uniprot database version
    with urllib.request.urlopen("https://unipept.ugent.be/private_api/metadata") as r:
        data = r.read()
    unipept_resp = json.loads(data)
    db_version = unipept_resp["db_version"]
    print(f"unipept database version: {db_version}", file=LOG_HANDLE)

    # map peptides
    trypticPeptides = []
    with open(peptide_file) as infile:
        for line in infile.readlines():
            content = line.strip()
            if content:
                trypticPeptides.append(content)

    # from unipept doku:
    # When performing bulk searches, we suggest splitting the input set over requests of 100 peptides each.
    chunksize = 100
    pep2prot = {}
    for chunk_start in range(0, len(trypticPeptides), chunksize):
        pep_chunk = trypticPeptides[chunk_start : chunk_start + chunksize]
        pep2prot.update(get_pep2uniprot_dict(pep_chunk))

    with open(out_path, "w") as outfile:
        json.dump(pep2prot, outfile)


def get_pep2uniprot_dict(pep_chunk):
    if len(pep_chunk) > 100:
        msg = f"{len(pep_chunk)} peptides provided. The maximum number of paptides to map whould not exceed 100. Performing request anyway."
        print(msg, file=LOG_HANDLE)
        raise Warning(msg)
    url = "http://api.unipept.ugent.be/api/v1/pept2prot.json"
    pep2prot = defaultdict(list)
    post_data = [("equate_il", "true")]
    post_data += [("input[]", x) for x in pep_chunk]
    data = urllib.parse.urlencode(post_data).encode()
    req = urllib.request.Request(url, data=data, method="POST")
    with urllib.request.urlopen(req) as response:
        rdata = response.read()
    unipept_resp = json.loads(rdata)
    for res in unipept_resp:
        pep2prot[res["peptide"]].append(res["uniprot_id"])
    return pep2prot


if "snakemake" in globals():
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        map_peptides_to_uniprot_ids(snakemake.input[0], snakemake.output[0])
