import json
import sys
import pandas as pd

LOG_HANDLE = sys.stderr


def map_taxids_to_uniprot_accessions(idmapping_selected_file: str, tax_ids_file: str, output_tax2acc_map: str,
                                     output_tax_ids: str):
    tax_ids = []
    with open(tax_ids_file) as infile:
        for line in infile.readlines():
            content = line.strip()
            if content:
                tax_ids.append(int(content))

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
    with open(output_tax2acc_map, "w") as outfile:
        json.dump(tax2acc_map, outfile, indent=4)
        outfile.write("\n")
    with open(output_tax_ids, "w") as outfile:
        outfile.write("\n".join(map(str, tax2acc_map.keys())))
        outfile.write("\n")


if snakemake := globals().get("snakemake"):
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        map_taxids_to_uniprot_accessions(
            idmapping_selected_file=snakemake.input["idmap"],
            tax_ids_file=snakemake.input["tax_ids"],
            output_tax2acc_map=snakemake.output["tax2acc_map"],
            output_tax_ids=snakemake.output["tax_ids"],
        )
