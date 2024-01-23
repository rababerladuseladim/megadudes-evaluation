from typing import Iterable

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["figure.figsize"] = [16 * 0.6, 9 * 0.6]
# context keywords: notebook, talk, paper, poster
sns.set_theme(context="talk", rc={"axes.grid": False})


def qc_plots(diamond_hits_file, input_fasta_file, output):
    df_hits = pd.read_table(
        diamond_hits_file,
        usecols=[0, 4],
        header=None,
        dtype={"query": str, "evalue": float},
        names=["query", "evalue"],
    )
    df_hits["length"] = df_hits["query"].str.len()

    with open(input_fasta_file) as fh:
        input_peptide_lengths = get_length_of_unique_sequences_from_fasta(fh)

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


def get_length_of_unique_sequences_from_fasta(fasta_handle: Iterable[str]) -> list[int]:
    sequences = []
    incomplete_sequence = ""
    for line in fasta_handle:
        if not line.strip():
            continue
        if line.startswith(">"):
            if incomplete_sequence and incomplete_sequence not in sequences:
                sequences.append(incomplete_sequence)
            incomplete_sequence = ""
            continue
        incomplete_sequence += line.strip()
    sequence_lengths = [len(seq) for seq in sequences]
    return sequence_lengths


def test_get_input_sequence_lengths():
    fasta = [
        ">foo",
        "ATGCCGCCG",
        ">bar",
        "ATGCCGCCG",
        "ATG",
        ">baz",
        "ATGCCGCCG "
        "",
    ]
    assert get_length_of_unique_sequences_from_fasta(fasta) == [9, 12]


def test_qc_plots(tmpdir):
    input_fasta_f = "/home/hennings/Projects/megadudes-evaluation/hpi-mnt/results/fastas/simulated_peptides_0.fasta"
    diamond_hits_f = "/home/hennings/Projects/megadudes-evaluation/hpi-mnt/results/diamond/simulated_peptides_0.tsv"
    output_plot = tmpdir / "hist_of_hits.svg"
    f = qc_plots(
        diamond_hits_file=diamond_hits_f,
        input_fasta_file=input_fasta_f,
        output=output_plot,
    )
    f.show()


if snakemake := globals().get("snakemake"):
    qc_plots(
        diamond_hits_file=snakemake.input[0],
        input_fasta_file=snakemake.input[1],
        output=snakemake.output[0],
    )
