import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["figure.figsize"] = [16 * 0.6, 9 * 0.6]
# context keywords: notebook, talk, paper, poster
sns.set_theme(context="talk", rc={"axes.grid": False})


def qc_plots(diamond_hits_file, input_fasta_file, output):
    df_hits = pd.read_table(
        diamond_hits_file,
        usecols=[0, 5, 7],
        header=None,
        dtype={"query": str, "identity": float, "evalue": float},
        names=["query", "identity", "evalue"],
    )
    df_hits["length"] = df_hits["query"].str.len()

    input_peptide_lengths = []
    with open(input_fasta_file) as fh:
        for line in fh:
            if line.startswith(">") or len(line.strip()) == 0:
                continue
            else:
                input_peptide_lengths.append(len(line.strip()))

    df_lengths = pd.concat(
        [
            pd.DataFrame({"length": df_hits.drop_duplicates(subset="query").loc[:, "length"],
                          "source": "hits"}),
            pd.DataFrame({"length": input_peptide_lengths,
                          "source": "input"}),
        ],
        ignore_index=True
    )

    # plot
    f, (ax0, ax1) = plt.subplots(2, 1, squeeze=True)
    sns.histplot(df_lengths, x="length", hue="source", ax=ax0)
    sns.histplot(df_hits, x="evalue", ax=ax1, log_scale=True)
    f.tight_layout()
    f.savefig(output)
    return f


def test_qc_plots(tmpdir):
    input_fasta_f = "/home/hennings/Projects/megadudes-evaluation/results/diamond/proc/peptides.fasta"
    diamond_hits_f = "/home/hennings/Projects/megadudes-evaluation/results/diamond/out.tsv"
    output_plot = tmpdir / "hist_of_hits.svg"
    f = qc_plots(
        diamond_hits_file=diamond_hits_f,
        input_fasta_file=input_fasta_f,
        output=output_plot,
    )
    f.show()


if "snakemake" in globals():
    qc_plots(
        diamond_hits_file=snakemake.input[0],
        input_fasta_file=snakemake.input[1],
        output=snakemake.output[0],
    )
