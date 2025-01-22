import json
import sys
import pandas as pd

LOG_HANDLE = sys.stderr


def map_taxids_to_uniprot_accessions(idmapping_selected_file: str, lineage_file: str, output_json: str):
    df_tax_ids = pd.read_csv(lineage_file, usecols=["query"], dtype={"query": int}, sep="\t")
    tax_ids = df_tax_ids["query"].to_list()

    df = pd.read_csv(
        idmapping_selected_file,
        compression="gzip" if idmapping_selected_file.endswith(".gz") else None,
        sep="\t",
        header=None,
        usecols=[0, 12],
        names=["acc", "taxid"],
        converters={12: lambda x: int(x) if int(x) in tax_ids else -1}
    )
    acc2taxid = df[df["taxid"] != -1]

    tax2acc_map = {}
    missing_tax_ids = []
    for tax_id in tax_ids:
        if tax_id in acc2taxid["taxid"].values:
            tax2acc_map[tax_id] = acc2taxid[acc2taxid["taxid"] == tax_id]["acc"].to_list()
        else:
            missing_tax_ids.append(tax_id)
    if missing_tax_ids:
        print(f"Taxids missing in idmapping file: {sorted(missing_tax_ids)}", file=LOG_HANDLE)
    with open(output_json, "w") as outfile:
        json.dump(tax2acc_map, outfile, indent=4)
        outfile.write("\n")


if snakemake := globals().get("snakemake"):
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        map_taxids_to_uniprot_accessions(
            idmapping_selected_file=snakemake.input["idmap"],
            lineage_file=snakemake.input["lineage"],
            output_json=snakemake.output[0],
        )
