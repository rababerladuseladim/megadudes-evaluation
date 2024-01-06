rule run_unipept_on_simulated_sample:
    input:
        "results/simulation/{sample}.txt",
    log:
        "logs/unipept/{sample}.txt",
    output:
        "results/unipept/{sample}.csv",
    conda:
        "../envs/ruby.yaml"
    shell:
        """
        gem install unipept
        cat {input} | prot2pept | peptfilter | unipept pept2lca -a -e > {output}
        """
