from typing import Iterable

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["figure.figsize"] = [16 * 0.6, 9 * 0.6]
# context keywords: notebook, talk, paper, poster
sns.set_theme(context="talk", rc={"axes.grid": False})


def qc_plots(diamond_hits_file, mmseqs2_hits_file, input_fasta_file, output):
    df_diamond_hits = pd.read_table(
        diamond_hits_file,
        usecols=[0, 4],
        header=None,
        dtype={"query": str, "evalue": float},
        names=["query", "evalue"],
    )
    df_diamond_hits["length"] = df_diamond_hits["query"].str.len()

    df_mmseqs2_hits = pd.read_table(
        mmseqs2_hits_file,
        usecols=[0, 4],
        header=None,
        dtype={"query": str, "evalue": float},
        names=["query", "evalue"],
    )
    df_mmseqs2_hits["length"] = df_mmseqs2_hits["query"].str.len()

    with open(input_fasta_file) as fh:
        input_peptide_lengths = get_length_of_unique_sequences_from_fasta(fh)

    df_lengths = pd.concat(
        [
            pd.DataFrame({"length": df_diamond_hits.drop_duplicates(subset="query").loc[:, "length"],
                          "source": "diamond"}),
            pd.DataFrame({"length": df_mmseqs2_hits.drop_duplicates(subset="query").loc[:, "length"],
                          "source": "mmseqs2"}),
            pd.DataFrame({"length": input_peptide_lengths,
                          "source": "input"}),
        ],
        ignore_index=True
    )

    df_evalue = pd.concat(
        [
            pd.DataFrame({"evalue": df_diamond_hits.loc[:, "evalue"],
                          "source": "diamond"}),
            pd.DataFrame({"evalue": df_mmseqs2_hits.loc[:, "evalue"],
                          "source": "mmseqs2"}),
        ],
        ignore_index=True
    )

    # plot
    f, (ax0, ax1) = plt.subplots(2, 1, squeeze=True)
    bins = list(range(df_lengths["length"].min(), df_lengths["length"].max()+1))
    sns.histplot(df_lengths, x="length", hue="source", ax=ax0, element="step", bins=bins)
    sns.histplot(df_evalue, x="evalue", hue="source", ax=ax1, log_scale=True, element="step")
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
    input_fasta_f = "/home/hennings/Projects/megadudes-evaluation/hpi-mnt/results/fastas/sample_kleiner_vs_kleiner.fasta"
    diamond_hits_f = "/home/hennings/Projects/megadudes-evaluation/hpi-mnt/results/diamond/sample_kleiner_vs_kleiner.tsv"
    mmseqs2_hits_f = "/home/hennings/Projects/megadudes-evaluation/hpi-mnt/results/mmseqs2/sample_kleiner_vs_kleiner.tsv"
    output_plot = tmpdir / "hist_of_hits.svg"
    f = qc_plots(
        diamond_hits_file=diamond_hits_f,
        mmseqs2_hits_file=mmseqs2_hits_f,
        input_fasta_file=input_fasta_f,
        output=output_plot,
    )
    f.show()


if snakemake := globals().get("snakemake"):
    qc_plots(
        diamond_hits_file=snakemake.input["diamond"],
        mmseqs2_hits_file=snakemake.input["mmseqs2"],
        input_fasta_file=snakemake.input["input"],
        output=snakemake.output[0],
    )
