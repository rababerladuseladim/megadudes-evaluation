import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["figure.figsize"] = [16 * 0.6, 9 * 0.6]
# context keywords: notebook, talk, paper, poster
sns.set_theme(context="talk", rc={"axes.grid": False})


def qc_plots(diamond_hits_file, output):
    df_hits = pd.read_table(
        diamond_hits_file,
        usecols=[0, 5],
        header=0,
        names=["query", "identity"],
    )
    df_hits["length"] = df_hits["query"].str.len()
    # plot
    f, axs = plt.subplots(1, 2, squeeze=False)
    sns.histplot(df_hits, x="length", ax=axs[0][0])
    sns.histplot(df_hits, x="identity", ax=axs[0][1])
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
