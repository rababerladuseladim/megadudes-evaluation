import pandas as pd


def map_uniprot_accessions_to_taxids(idmapping_selected_file, accession_file, output_file, log_handle):
    accessions = []
    with open(accession_file) as infile:
        for line in infile.readlines():
            content = line.strip()
            if content:
                accessions.append(content)
    accessions = set(accessions)

    df = pd.read_csv(
        idmapping_selected_file,
        compression="gzip" if idmapping_selected_file.endswith(".gz") else None,
        sep="\t",
        header=None,
        usecols=[0, 12],
        names=["acc", "taxid"],
        converters={0: lambda x: x if x in accessions else "missing"}
    )
    acc2taxid = df[df["acc"] != "missing"]

    taxids = []
    missing_accessions = []
    for acc in accessions:
        if acc in acc2taxid["acc"].values:
            taxids += acc2taxid[acc2taxid["acc"] == acc]["taxid"].to_list()
        else:
            missing_accessions.append(acc)
    if missing_accessions:
        print(f"Accessions missing in idmapping file: {missing_accessions}", file=log_handle)
    taxids = sorted(set(taxids))
    with open(output_file, "w") as outfile:
        outfile.write("\n".join(map(str, taxids)))
        outfile.write("\n")


if snakemake := globals().get("snakemake"):
    with open(snakemake.log[0], "w") as log_handle:
        map_uniprot_accessions_to_taxids(
            idmapping_selected_file=snakemake.input["idmap"],
            accession_file=snakemake.input["accessions"],
            output_file=snakemake.output[0],
            log_handle=log_handle
        )
