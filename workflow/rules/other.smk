# An example collection of Snakemake rules imported in the main Snakefile.
rule map_peptides_to_uniprot_ids:
    input:
        config["peptides"]
    output:
        "results/map_peptides_to_uniprot_ids/pep2uniprot_ids.json"
    run:
        from collections import defaultdict
        import json
        import urllib.parse
        import urllib.request
        url = 'http://api.unipept.ugent.be/api/v1/pept2prot.json'
        trypticPeptides = []
        with open(input[0]) as infile:
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

        with open(output[0],"w") as outfile:
            json.dump(pep2prot, outfile)
