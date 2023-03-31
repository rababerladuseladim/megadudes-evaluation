rule plot_megadudes_qc:
    input:
        ground_truth=config["ground_truth"],
        unipept_results=[config["unipept_result"]],
        megadudes_results=[
            "results/megadudes/npz.out",
            "results/megadudes/blast.out",
        ],
    log:
        "logs/megadudes/plot_megadudes_qc.txt",
    output:
        report("plots/megadudes/qc.svg"),
    conda:
        "../envs/qc_plots.yaml"
    script:
        "../scripts/plot_megadudes_qc.py"