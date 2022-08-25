import os.path

rule convert_gzip_to_block_gzip:
    input: "{f}.gz"
    output: "{f}.bgz"
    conda: "../envs/samtools.yaml"
    log:
        zcat="logs/megadudes/convert_gzip_to_block_gzip-zcat-{f}.txt",
        bgzip="logs/megadudes/convert_gzip_to_block_gzip-bgzip-{f}.txt",
    shell: "zcat {input} 2>{log.zcat}| bgzip -c 2>{log.bgzip}> {output}"


rule map_peptides_from_uniparc_to_uniprot_ids:
    input:
        pep_tsv=config["msfragger_peptides_tsv"],
        dudes_db="results/megadudes/proc/dudes_db.npz",
        uniparc_fasta="resources/uniprot/uniparc.fasta.bgz",
        idmap="resources/uniprot/idmapping_selected.tab.gz"
    log:
        "logs/megadudes/map_peptides_from_uniparc_to_uniprot_ids.txt"
    output:
        "results/megadudes/proc/mapped_peptides.npz"
    conda:
        "../envs/megadudes.yaml"
    shell:
        """
        dudes-parse-peptides \
        -m {input.pep_tsv} \
        -d {input.dudes_db} \
        -f {input.uniparc_fasta} \
        -i {input.idmap} \
        -o {output} \
        >> {log}
        """



rule build_megadudes_db:
    input:
        idmap="resources/uniprot/idmapping_selected.tab.gz",
        uniprot_fastas=expand("resources/uniprot/{db}.fasta.gz",
            db=[
            "swissprot",
            "trembl"
        ]),
        ncbi_nodes="resources/ncbi/nodes.dmp",
        ncbi_names="resources/ncbi/names.dmp"
    log:
        "logs/megadudes/build_megadudes_db.txt"
    output:
        dudes_db="results/megadudes/proc/dudes_db.npz"
    conda:
        "../envs/megadudes.yaml"
    threads: 99
    params:
        db_base_name=lambda wildcards, output: os.path.splitext(output["dudes_db"])[0]
    shell:
        """
        dudesdb -m up \
        -f {input.uniprot_fastas} \
        -g {input.idmap} \
        -n {input.ncbi_nodes} \
        -a {input.ncbi_names} \
        -o {params.db_base_name} \
        -t {threads} \
        >> {log}
        """


rule run_megadudes:
    input:
        dudes_db="results/megadudes/proc/dudes_db.npz",
        pep_db="results/megadudes/proc/mapped_peptides.npz"
    log:
        "logs/megadudes/run_megadudes-normalize_{normalize}.txt"
    output:
        table="results/megadudes/result-normalize_{normalize}.out",
        plots=directory("plots/megadudes-normalize_{normalize}")
    wildcard_constraints:
        normalize="true|false"
    conda:
        "../envs/megadudes.yaml"
    threads:
        99
    params:
        out_wo_ext=lambda wildcards, output: os.path.splitext(output.table)[0],
        additional_args=lambda w: [] if w.normalize == "true" else ["--no-normalize"]
    shell:
        """
        dudes \
        -s {input.pep_db} \
        -d {input.dudes_db} \
        -t {threads} \
        -o {params.out_wo_ext}\
        {params.additional_args} \
        --debug_plots_dir {output.plots} \
        --debug 2> {log}
        """


rule run_megadudes_on_diamond_sam:
    input:
        dudes_db="results/megadudes/proc/dudes_db.npz",
        sam_file="results/diamond/out.sam"
    log:
        "logs/megadudes/run_megadudes-normalize_{normalize}-diamond.txt"
    output:
        table="results/megadudes/result-normalize_{normalize}-diamond.out",
        plots=directory("plots/megadudes-normalize_{normalize}-diamond")
    wildcard_constraints:
        normalize="true|false"
    conda:
        "../envs/megadudes.yaml"
    threads:
        99
    params:
        out_wo_ext=lambda wildcards, output: os.path.splitext(output.table)[0],
        additional_args=lambda w: [] if w.normalize == "true" else ["--no-normalize"]
    shell:
        """
        dudes \
        -s {input.sam_file} \
        -d {input.dudes_db} \
        -t {threads} \
        -o {params.out_wo_ext}\
        {params.additional_args} \
        --debug_plots_dir {output.plots} \
        --debug 2> {log}
        """
