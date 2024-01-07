rule install_unipept:
    log:
        "logs/unipept/install.txt",
    output:
        "results/unipept/install_success",
    conda:
        "../envs/ruby.yaml"
    shell:
        """
        gem install unipept > {log} 2>&1;
        ln -s $(which ruby) $(dirname $(which unipept))/ruby;  # unipept looks for ruby in the gem bin dir
        touch {output};
        """


rule run_unipept_on_simulated_sample:
    input:
        "results/unipept/install_success",
        "results/simulation/{sample}.txt",
    log:
        "logs/unipept/{sample}.txt",
    output:
        "results/unipept/{sample}.csv",
    conda:
        "../envs/ruby.yaml"
    shell:
        """
        cat {input} | prot2pept | peptfilter | unipept pept2lca -a -e > {output} 2>{log};
        """
