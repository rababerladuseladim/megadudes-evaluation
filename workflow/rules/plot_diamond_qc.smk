rule plot_diamond_qc:
    input:
        "results/diamond/out.tsv",
    output:
        report("plots/diamond/qc.svg"),
    log:
        "logs/diamond/plot_diamond_qc.txt",
    conda:
        "../envs/qc_plots.yaml"
    script:
        "../scripts/plot_diamond_qc.py"
