rule qc_plots:
    input:
        ground_truth="resources/kleiner_ground_truth-equal_protein.csv",
        unipept_results=[],
        megadudes_results=["results/megadudes/result.out"]
    log:
        "logs/qc_plots.txt"
    output:
        plots="plots/qc.svg"
    conda:
        "../envs/qc_plots.yaml"
    script:
        "../scripts/qc_plots.py"
