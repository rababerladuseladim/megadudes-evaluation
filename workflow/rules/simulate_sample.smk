rule sample_taxons:
    input:
        "resources/uniprot/speclist.txt",
    output:
        "results/simulation/sample_taxons.txt",
    log:
        "logs/simulation/sample_taxons.txt",
    script:
        "../scripts/sample_taxons.py"


rule map_taxids_to_uniprot_accessions:
    input:
        idmap="resources/uniprot/idmapping_selected.tab.gz",
        taxids="results/simulation/sample_taxons.txt",
    output:
        "results/simulation/tax2accessions.json",
    log:
        "logs/simulation/map_taxids_to_uniprot_accssions.txt",
    script:
        "../scripts/map_taxids_to_uniprot_accessions.py"


rule sample_peptides:
    input:
        accessions="results/simulation/tax2accessions.json",
    output:
        "results/simulation/peptides.txt",
    log:
        "logs/simulation/sample_peptides.txt",
    conda:
        "../envs/pyteomics.yaml"
    script:
        "../scripts/sample_peptides.py"


rule convert_simulated_peptides_to_fasta:
    input:
        "results/simulation/peptides_{N}.txt",
    output:
        "results/fastas/simulated_peptides_{N}.fasta",
    log:
        "logs/simulation/convert_peptides_txt_to_fasta_{N}.txt",
    shell:
        "cat {input} | sed 's/.*/>&\\n&/' > {output} 2>{log}"
