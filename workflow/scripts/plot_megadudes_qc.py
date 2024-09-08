from pathlib import Path

import pandas as pd

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake import script


try:
    import seaborn as sns
    import matplotlib.pyplot as plt
except ImportError:
    pass


TAX_LEVELS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies"]


def read_ground_truth_file(ground_truth_file):
    df_t_taxids = pd.read_csv(
        ground_truth_file,
        usecols=TAX_LEVELS,
        sep="\t",
        dtype={
            c: pd.Int64Dtype()  # support nan values
            for c in TAX_LEVELS
        },
    )
    return df_t_taxids


def get_value_overlap(hits: pd.Series, ground_truth: pd.Series):
    """
    build series of values from both input series and annotate them with TP, FP, FN

    :param hits:
    :param ground_truth:
    :return: series with index: union of values from hits and ground_truth, values: TP if index in both,
        FP if only in hits, FN if only in ground_truth
    """
    left = {t for t in hits if not (pd.isna(t) or (t == "-"))}
    right = {t for t in ground_truth if not (pd.isna(t) or (t == "-"))}
    taxa = left.union(right)
    dct = dict()
    for t in taxa:
        found = t in left
        expected = t in right
        if found and expected:
            dct[t] = "TP"
        elif found:
            dct[t] = "FP"
        elif expected:
            dct[t] = "FN"
        else:
            raise KeyError(f"Taxon '{t}' missing in both Series!")
    return pd.Series(dct)


def get_diamond_hit_counts(f: Path | str, ground_truth_df_tax_ids: pd.DataFrame):
    df_dmnd = read_ground_truth_file(f)
    df = pd.DataFrame(columns=TAX_LEVELS, index=["TP", "FP", "FN"],)
    for t in TAX_LEVELS:
        s = get_value_overlap(df_dmnd[t], ground_truth_df_tax_ids[t])
        df[t] = s.value_counts()
    df = df.fillna(0)
    df["eval"] = df.index
    return df


def get_unipept_hit_counts(f, ground_truth_df_tax_ids):
    df_uni = pd.read_csv(
        f,
        usecols=[f"{t}_id" for t in TAX_LEVELS],
        dtype=pd.Int64Dtype()
    )
    df_uni.rename(lambda x: x.strip("_id"), axis="columns", inplace=True)
    df = pd.DataFrame(columns=TAX_LEVELS, index=["TP", "FP", "FN"],)
    for t in TAX_LEVELS:
        s = get_value_overlap(df_uni[t], ground_truth_df_tax_ids[t])
        df[t] = s.value_counts()
    df = df.fillna(0)
    df["eval"] = df.index
    return df


def get_megadudes_hit_counts(f, ground_truth_df_tax_ids):
    df_mega = pd.read_csv(f, sep="\t", skiprows=5)
    df_mega = df_mega.rename(lambda c: c.lower().strip("@"), axis=1)
    df = pd.DataFrame(columns=TAX_LEVELS, index=["TP", "FP", "FN"])
    for t in TAX_LEVELS:
        s_mega = (
            df_mega.loc[df_mega["rank"] == t, "taxpath"]
            .map(lambda x: x.split("|")[-1])
            .astype(int)
        )
        s = get_value_overlap(s_mega, ground_truth_df_tax_ids[t])
        s_classification = s.value_counts()
        for k in ["TP", "FP", "FN"]:
            if k not in s_classification:
                s_classification[k] = 0
        df[t] = s_classification
    df = df.fillna(0)
    df["eval"] = df.index
    return df


def calc_eval_metrics(df_hits: pd.DataFrame) -> pd.DataFrame:
    df_hits_long = df_hits.melt(id_vars=["method", "eval"], var_name="taxon_level")
    df_hits_long["value"] = df_hits_long["value"].apply(int)
    grouped = df_hits_long.groupby(["method", "taxon_level"], sort=False)

    df_sens_g = grouped.apply(
        lambda df: df[df["eval"] == "TP"]["value"].values[0]
        * 100
        / df[df["eval"].isin(["TP", "FN"])]["value"].sum()
    )
    # specificity: TN / (TN + FP)
    # how to access TN?
    # precision: TP / (TP + FP)
    df_precision_g = grouped.apply(
        lambda df: df[df["eval"] == "TP"]["value"].values[0]
        * 100
        / df[df["eval"].isin(["TP", "FP"])]["value"].sum()
    )
    df_f_one_g = (2 * df_precision_g * df_sens_g) / (df_precision_g + df_sens_g)
    df_fdr_g = grouped.apply(
        lambda df: df[df["eval"] == "FP"]["value"].values[0]
        * 100
        / df[df["eval"].isin(["TP", "FP"])]["value"].sum()
    )
    df_eval = pd.DataFrame(
        {
            "sensitivity": df_sens_g,
            "precision": df_precision_g,
            "F1-score": df_f_one_g,
            "fdr": df_fdr_g,
        }
    )
    df_eval = df_eval.reset_index()
    return df_eval


def plot_qc(df_plt, output):
    metrics = ["sensitivity", "precision", "F1-score"]
    f, axes = plt.subplots(ncols=3, figsize=(15, 9), sharex="all", sharey="all")
    for metric, ax in zip(metrics, axes):
        sns.lineplot(
            data=df_plt,
            x="taxon_level",
            y=metric,
            hue="method",
            ax=ax,
            # legend only in first plot
            legend=bool(metric == metrics[0]),
            marker="s",
        )
        ax.set_title(metric)
        ax.set_ylim(0, 105)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        ax.xaxis.label.set_visible(False)
    axes[0].yaxis.label.set_visible(False)
    axes[-1].tick_params(labelright=True)
    f.savefig(output)
    return f


def qc_plots(
    ground_truth_file: str,
    unipept_file: str,
    megadudes_file_list: list[str],
    output: str,
    diamond_file: str | None = None
):
    # context keywords: notebook, talk, paper, poster
    sns.set_theme(context="talk", rc={"axes.grid": False})
    plt.rcParams["figure.figsize"] = [16 * 0.6, 9 * 0.6]
    hits: list[pd.DataFrame] = []
    # read ground truth
    df_gt_taxids = read_ground_truth_file(ground_truth_file)
    # get diamond hits
    df_dmnd = get_diamond_hit_counts(diamond_file, df_gt_taxids) if diamond_file else pd.DataFrame()
    if not df_dmnd.empty:
        df_dmnd["method"] = "diamond"
    hits.append(df_dmnd)
    # get unipept hits
    df_uni = get_unipept_hit_counts(unipept_file, df_gt_taxids)
    df_uni["method"] = "unipept"
    hits.append(df_uni)
    # get megadudes hits
    for megadudes_file in megadudes_file_list:
        df_mega = get_megadudes_hit_counts(
                megadudes_file, df_gt_taxids
            )
        df_mega["method"] = "megadudes on " + Path(megadudes_file).parent.name
        hits.append(df_mega)
    # build TP/FP/FN dataframe
    df_hits = pd.concat(hits, ignore_index=True)
    # build evaluation dataframe, containing sensitivity, precision, f1-score, fdr
    df_eval = calc_eval_metrics(df_hits)
    # plot
    # drop subspecies, as no method has any true hits on this level
    df_plt = df_eval[df_eval["taxon_level"] != "subspecies"]
    return plot_qc(df_plt, output), df_eval


def test_qc_plots(tmpdir):
    ground_truth = "/home/hennings/Projects/megadudes-evaluation/resources/kleiner_ground_truth-equal_protein-full_taxid_lineage.csv"
    unipept_file = "/home/hennings/Nextcloud/Notes/project support material/taxfdr/files/01-researchgap/02-tax_analysis/unipept_results-Kleiner_P_msfragger.csv"
    megadudes_file = "/home/hennings/Nextcloud/Notes/project support material/taxfdr/files/02-megadudes/v0.9/megadudes-result.out"
    output_plot = tmpdir / "qc_plot.svg"
    print(output_plot)
    f, df_eval = qc_plots(
        ground_truth_file=ground_truth,
        unipept_file=unipept_file,
        megadudes_file_list=[megadudes_file],
        output=output_plot,
    )
    f.show()


snakemake: "script.Snakemake"
if snakemake := globals().get("snakemake"):
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        qc_plots(
            ground_truth_file=snakemake.input.ground_truth,
            diamond_file=snakemake.input.get("diamond_result"),
            unipept_file=snakemake.input.unipept_result,
            megadudes_file_list=snakemake.input.megadudes_results,
            output=snakemake.output[0],
        )
