rule taxon_identification_qc_plots:
    input:
        ground_truth=config["ground_truth"],
        unipept_results=[config["unipept_result"]],
        megadudes_results=[
            "results/megadudes/npz-normalize_false.out",
            "results/megadudes/blast-normalize_false.out",
        ],
    log:
        "logs/qc_plots.txt",
    output:
        report("plots/taxon_identification_qc.svg"),
    conda:
        "../envs/qc_plots.yaml"
    script:
        "../scripts/taxon_identification_qc_plots.py"
