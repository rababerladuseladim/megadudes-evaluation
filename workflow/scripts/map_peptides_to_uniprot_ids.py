from collections import defaultdict
import json
import urllib.parse
import urllib.request


def do_something(peptide_file, out_path, log_path):
    url = 'http://api.unipept.ugent.be/api/v1/pept2prot.json'
    trypticPeptides = []
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
    do_something(snakemake.input[0], snakemake.output[0], snakemake.log)
