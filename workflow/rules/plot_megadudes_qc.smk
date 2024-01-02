rule plot_megadudes_qc_sample:
    input:
        ground_truth=lambda wc: samples.loc[wc.sample_name, "ground_truth"],
        unipept_results=lambda wc: [samples.loc[wc.sample_name, "unipept_result"]],
        megadudes_results=["results/megadudes/sample_{sample_name}.out"],
    log:
        "logs/megadudes/plot_megadudes_qc-sample_{sample_name}.txt",
    output:
        report("plots/megadudes/qc-sample_{sample_name}.svg"),
    conda:
        "../envs/qc_plots.yaml"
    script:
        "../scripts/plot_megadudes_qc.py"


rule plot_megadudes_qc_simulation:
    input:
        ground_truth="results/simulation/sample_taxons_lineage_{repeat}.tsv",
        unipept_results=[],
        megadudes_results=["results/megadudes/simulated_peptides_{repeat}.out"],
    log:
        "logs/megadudes/plot_megadudes_qc-simulated_peptides_{repeat}.txt",
    output:
        report("plots/megadudes/qc-simulated_peptides_{repeat}.svg"),
    conda:
        "../envs/qc_plots.yaml"
    script:
        "../scripts/plot_megadudes_qc.py"
