import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["figure.figsize"] = [16 * 0.6, 9 * 0.6]
# context keywords: notebook, talk, paper, poster
sns.set_theme(context="talk", rc={"axes.grid": False})


def qc_plots(diamond_hits_file, output):
    df_hits = pd.read_table(
        diamond_hits_file,
        usecols=[0, 5, 7],
        header=0,
        names=["query", "identity", "evalue"],
    )
    df_hits["length"] = df_hits["query"].str.len()
    # plot
    f, axs = plt.subplots(3, 1, squeeze=False)
    sns.histplot(df_hits, x="length", ax=axs[0][0])
    sns.histplot(df_hits, x="identity", ax=axs[1][0])
    sns.histplot(df_hits, x="evalue", ax=axs[2][0], log_scale=True)
    f.tight_layout()
    f.savefig(output)
    return f


def test_qc_plots(tmpdir):
    diamond_hits_f = "/home/hennings/Projects/megadudes-evaluation/results/diamond/out.tsv"
    output_plot = tmpdir / "hist_of_hits.svg"
    f = qc_plots(
        diamond_hits_file=diamond_hits_f,
        output=output_plot,
    )
    f.show()


if "snakemake" in globals():
    qc_plots(
        diamond_hits_file=snakemake.input[0],
        output=snakemake.output[0],
    )
