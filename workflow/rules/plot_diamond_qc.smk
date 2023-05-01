rule plot_diamond_qc:
    input:
        "results/diamond/{sample}.tsv",
        "results/fastas/{sample}.fasta",
    output:
        report("plots/diamond/qc-{sample}.svg"),
    log:
        "logs/diamond/plot_diamond_qc-{sample}.txt",
    conda:
        "../envs/qc_plots.yaml"
    script:
        "../scripts/plot_diamond_qc.py"
