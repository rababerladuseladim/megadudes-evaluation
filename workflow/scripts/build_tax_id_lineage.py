import pandas as pd
from common import TAX_LEVELS
from pathlib import Path


def parse_nodes(nodes_path: str) -> pd.DataFrame:
    """Read ncbi nodes.dmp file into pandas DataFrame.

    Args:
        nodes_path: Path to ncbi nodes.dmp file

    Returns:
        DataFrame with the columns 'tax_id', 'parent_tax_id', and 'rank'
    """
    return pd.read_table(
        nodes_path,
        sep=r"\t\|\t",
        engine="python",
        names=["tax_id", "parent_tax_id", "rank"],
        index_col=0,
        usecols=[0, 1, 2],
    )


def build_tax_id_lineage(tax_id: int, nodes: pd.DataFrame) -> dict[str, int]:
    """Build lineage of tax_id and all parent tax_ids.

    lineage is restricted to tax_id that have a rank which is in TAX_LEVELS

    Args:
        tax_id: tax_id
        nodes: nodes dataframe, expected to have tax_id as index, first columns are parent tax_id and secondo columns
            are rank

    Returns:
        lineage as mapping from rank to tax_id
    """
    lineage: dict[str, int] = {}
    while True:
        parent_tax_id, rank = nodes.loc[tax_id, :].values
        lineage[rank] = tax_id
        if tax_id == parent_tax_id:
            break
        tax_id = parent_tax_id
    return {rank: tax_id for rank, tax_id in lineage.items() if rank in TAX_LEVELS}


def build_tax_id_lineage_tsv(tax_id_files: list[str], nodes: str) -> None:
    """Create one tsv file containing TAX_LEVEL lineage for each tax_id_file."""
    nodes = parse_nodes(nodes)
    for file_ in map(Path, tax_id_files):
        with file_.open() as f:
            tax_ids = [int(l) for l in f]
        lineage = [build_tax_id_lineage(i, nodes) for i in tax_ids]
        df = pd.DataFrame(lineage, columns=TAX_LEVELS)
        f_out = file_.parent / (file_.stem + "_lineage.tsv")
        df.to_csv(f_out, sep="\t", index=False)


def test_build_tax_id_lineage_tsv():
    simulated_taxon_file = "/home/hennings/Projects/megadudes-evaluation/test/unit/simulate_sample/build_tax_id_lineage/data/results/simulation/sample_taxons_1.txt"
    nodes_path = "/home/hennings/Projects/megadudes-evaluation/test/unit/simulate_sample/build_tax_id_lineage/data/resources/ncbi/nodes.dmp"
    build_tax_id_lineage_tsv([simulated_taxon_file], nodes_path)


def test_create_nodes_dmp() -> None:
    """Create minimal nodes.dmp with only the taxids from the lineageas in the provided taxon file"""
    simulated_taxon_file = "/home/hennings/Projects/megadudes-evaluation/test/unit/simulate_sample/sample_taxons/expected/results/simulation/sample_taxons_1.txt"
    nodes_path = "/home/hennings/Projects/megadudes-evaluation/resources/ncbi/nodes.dmp"
    test_db = "/home/hennings/Projects/megadudes-evaluation/test/unit/simulate_sample/build_tax_id_lineage/data/resources/ncbi/nodes.dmp"
    nodes = pd.read_table(
        nodes_path,
        sep=r"\t\|\t",
        header=None,
        engine="python",
        index_col=0,
    )
    nodes.drop(nodes.columns[-1], axis=1, inplace=True)

    with open(simulated_taxon_file) as f:
        tax_ids = [int(l) for l in f]
    lines = []
    for tax_id in tax_ids:
        while True:
            line = nodes.loc[tax_id, :]
            parent_tax_id = line.iloc[0]
            lines.append(line)
            if tax_id == parent_tax_id:
                break
            tax_id = parent_tax_id
    with (open(test_db, "w") as f):
        last_name = ""
        for line in sorted(lines, key=lambda l: l.name):
            if line.name == last_name:
                continue
            f.write("\t|\t".join(map(str, [line.name, *line.to_list()])) + "\t|\n")
            last_name = line.name


if "snakemake" in globals():
    build_tax_id_lineage_tsv(tax_id_files=snakemake.input["tax_id_files"], nodes=snakemake.input["ncbi_nodes"])
