rule sample_taxons:
    input:
        "resources/uniprot/speclist.txt",
    output:
        "results/sample_taxons/sample_taxons.txt",
    log:
        "logs/sample_taxons.txt",
    script:
        "../scripts/sample_taxons.py"
