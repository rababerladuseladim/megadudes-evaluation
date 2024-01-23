rule plot_megadudes_qc_sample:
    input:
        ground_truth=lambda wc: samples.loc[wc.sample_name, "ground_truth"],
        unipept_result=lambda wc: samples.loc[wc.sample_name, "unipept_result"],
        megadudes_result="results/megadudes/sample_{sample_name}.out",
    log:
        "logs/megadudes/plot_megadudes_qc-sample_{sample_name}.txt",
    output:
        report(
            "plots/megadudes/qc-sample_{sample_name}.svg",
            category="qc",
            subcategory="megadudes",
        ),
    conda:
        "../envs/qc_plots.yaml"
    script:
        "../scripts/plot_megadudes_qc.py"


rule plot_megadudes_qc_simulation:
    input:
        ground_truth="results/simulation/sample_taxons_lineage_{repeat}.tsv",
        diamond_result="results/diamond/simulated_peptides_{repeat}-lineage.tsv",
        unipept_result="results/unipept/peptides_{repeat}.csv",
        megadudes_result="results/megadudes/simulated_peptides_{repeat}.out",
    log:
        "logs/megadudes/plot_megadudes_qc-simulated_peptides_{repeat}.txt",
    output:
        report(
            "plots/megadudes/qc-simulated_peptides_{repeat}.svg",
            category="qc",
            subcategory="megadudes",
        ),
    conda:
        "../envs/qc_plots.yaml"
    script:
        "../scripts/plot_megadudes_qc.py"
