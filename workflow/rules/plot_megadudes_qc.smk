rule plot_megadudes_qc:
    input:
        ground_truth=config["ground_truth"],
        unipept_results=[config["unipept_result"]],
        megadudes_results=[
            "results/megadudes/result.out",
        ],
    output:
        report("plots/megadudes/qc.svg"),
    log:
        "logs/megadudes/plot_megadudes_qc.txt",
    conda:
        "../envs/qc_plots.yaml"
    script:
        "../scripts/plot_megadudes_qc.py"
