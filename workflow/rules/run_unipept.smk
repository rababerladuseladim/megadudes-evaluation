rule install_unipept:
    log:
        "logs/unipept/install.txt",
    output:
        "results/unipept/install_success",
    conda:
        "../envs/ruby.yaml"
    shell:
        """
        gem install unipept -v 3.1.0 > {log} 2>&1;
        ln -s $(which ruby) $(dirname $(which unipept))/ruby;  # unipept looks for ruby in the gem bin dir
        touch {output};
        """


rule run_unipept:
    input:
        install_success="results/unipept/install_success",
        peptides="results/peptides/{sample}.txt",
    log:
        "logs/unipept/{sample}.txt",
    output:
        "results/unipept/{sample}.csv",
    conda:
        "../envs/ruby.yaml"
    shell:
        """
        cat {input.peptides} | prot2pept | peptfilter | unipept pept2lca -a -e > {output} 2>{log};
        """
