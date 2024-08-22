import os.path


rule build_megadudes_db:
    input:
        idmap="resources/uniprot/idmapping_selected.tab.gz",
        uniprot_fastas=config["query_dbs"],
        ncbi_nodes="resources/ncbi/nodes.dmp",
        ncbi_names="resources/ncbi/names.dmp",
    log:
        stdout="logs/megadudes/build_megadudes_db-stdout.txt",
        stderr="logs/megadudes/build_megadudes_db-stderr.txt",
    output:
        dudes_db="results/megadudes/proc/dudes_db.npz",
    conda:
        "../envs/megadudes.yaml"
    threads: 99
    params:
        db_base_name=lambda wildcards, output: os.path.splitext(output["dudes_db"])[0],
    shell:
        """
        dudesdb -m up \
        -f {input.uniprot_fastas} \
        -g {input.idmap} \
        -n {input.ncbi_nodes} \
        -a {input.ncbi_names} \
        -o {params.db_base_name} \
        -t {threads} > {log.stdout} 2> {log.stderr}
        """


rule run_megadudes:
    input:
        dudes_db="results/megadudes/proc/dudes_db.npz",
        mmseqs2_result="results/mmseqs2_top_10/{sample}.tsv",
    log:
        stdout="logs/megadudes/run_megadudes-{sample}-stdout.txt",
        stderr="logs/megadudes/run_megadudes-{sample}-stderr.txt",
    output:
        result=report("results/megadudes/{sample}.out", category="megadudes"),
        plots=directory("plots/megadudes/{sample}/scores"),
    conda:
        "../envs/megadudes.yaml"
    threads: 99
    params:
        result_wo_ext=lambda wildcards, output: os.path.splitext(output.result)[0],
    shell:
        """
        dudes \
        -c {input.mmseqs2_result} \
        -d {input.dudes_db} \
        -t {threads} \
        -o {params.result_wo_ext}\
        --debug_plots_dir {output.plots} \
        --debug > {log.stdout} 2> {log.stderr}
        """
