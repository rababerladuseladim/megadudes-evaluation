import pandas as pd
from common import TAX_LEVELS
import sys

LOG_HANDLE = sys.stderr


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
        nodes: nodes dataframe, expected to have tax_id as index, first column contains the parent tax_id and the second
            column the rank

    Returns:
        lineage as mapping from rank to tax_id
    """
    lineage: dict[str, int] = {}
    while True:
        parent_tax_id, rank = nodes.loc[tax_id, :].values
        if rank in TAX_LEVELS:
            lineage[rank] = tax_id
        if tax_id == parent_tax_id:
            break
        tax_id = parent_tax_id
    return lineage


def filter_tax_ids_and_build_lineage_tsv(
    tax_ids_path: str, nodes_path: str, output_lineage_tsv: str
) -> None:
    """Filter tax_ids and build lineage based on ncbi nodes.dmp.

    Filter out tax_ids not in nodes.dmp, write remaining to a new file.
    Create a lineage table containing the tax_ids parents. Parents are restricted to ranks present in TAX_LEVELS.

    Args:
        tax_ids_path: path to file containing one tax_id per line
        nodes_path: path to ncbi nodes.dmp file
        output_lineage_tsv: path to output lineage tsv file
    """
    nodes = parse_nodes(nodes_path)
    with open(tax_ids_path) as f:
        tax_ids = [int(l) for l in f]
    lineage: list[dict[str, int]] = []
    for tax_id in tax_ids:
        try:
            lineage_current_tax_id = build_tax_id_lineage(tax_id, nodes)
        except KeyError:
            print(f"Dropping Taxonomy ID {tax_id} because it is not in nodes.dmp", file=LOG_HANDLE)
            continue
        lineage.append({**lineage_current_tax_id, "query": tax_id})

    df = pd.DataFrame(lineage, columns=[*TAX_LEVELS, "query"])
    df.to_csv(output_lineage_tsv, sep="\t", index=False)


def test_build_tax_id_lineage_tsv():
    simulated_taxon_file = "/home/hennings/Projects/megadudes-evaluation/test/rules/simulate_sample/build_tax_id_lineage/data/results/simulation/sample_taxons_1.txt"
    nodes_path = "/home/hennings/Projects/megadudes-evaluation/test/rules/simulate_sample/build_tax_id_lineage/data/resources/ncbi/nodes.dmp"
    filter_tax_ids_and_build_lineage_tsv(simulated_taxon_file, nodes_path)


def test_create_nodes_dmp() -> None:
    """Create minimal nodes.dmp with only the taxids from the lineage as in the provided taxon file"""
    simulated_taxon_file = "/test/unit/simulate_sample/sample_taxons/expected/results/simulation/sample_taxons_lineage_1.txt"
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
    with open(test_db, "w") as f:
        last_name = ""
        for line in sorted(lines, key=lambda l: l.name):
            if line.name == last_name:
                continue
            f.write("\t|\t".join(map(str, [line.name, *line.to_list()])) + "\t|\n")
            last_name = line.name


if snakemake := globals().get("snakemake"):
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        filter_tax_ids_and_build_lineage_tsv(
            tax_ids_path=snakemake.input["tax_ids"],
            nodes_path=snakemake.input["ncbi_nodes"],
            output_lineage_tsv=snakemake.output[0],
        )
