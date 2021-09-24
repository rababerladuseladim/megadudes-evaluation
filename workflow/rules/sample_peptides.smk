rule sample_peptides:
    input:
        accessions="results/accessions.txt"
    output:
        "results/sample_peptides/peptides.txt"
    log:
        "logs/sample_peptides.txt"
    script:
        "../scripts/sample_peptides.py"
