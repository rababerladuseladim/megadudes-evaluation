import json

import pandas as pd


def map_taxids_to_uniprot_accessions(idmapping_selected_file, taxid_file, output_json, log_handle):
    taxids = []
    with open(taxid_file) as infile:
        for line in infile.readlines():
            content = line.strip()
            if content:
                taxids.append(int(content))

    df = pd.read_csv(
        idmapping_selected_file,
        compression="gzip" if idmapping_selected_file.endswith(".gz") else None,
        sep="\t",
        header=None,
        usecols=[0, 12],
        names=["acc", "taxid"],
        converters={12: lambda x: int(x) if int(x) in taxids else -1}
    )
    acc2taxid = df[df["taxid"] != -1]

    tax2acc_map = {}
    missing_taxids = []
    for taxid in taxids:
        if taxid in acc2taxid["taxid"].values:
            tax2acc_map[taxid] = acc2taxid[acc2taxid["taxid"] == taxid]["acc"].to_list()
        else:
            missing_taxids.append(taxid)
    if missing_taxids:
        print(f"Taxids missing in idmapping file: {missing_taxids}", file=log_handle)
    with open(output_json, "w") as outfile:
        json.dump(tax2acc_map, outfile)
        outfile.write("\n")


def main():
    with open(snakemake.log[0], "w") as log_handle:
        map_taxids_to_uniprot_accessions(
            idmapping_selected_file=snakemake.input["idmap"],
            taxid_file=snakemake.input["taxids"],
            output_json=snakemake.output[0],
            log_handle=log_handle
        )


if "snakemake" in globals():
    main()
