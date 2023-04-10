rule sample_taxons:
    input:
        "resources/uniprot/speclist.txt",
    log:
        "logs/sample_taxons.txt",
    output:
        directory("results/sample_taxons/"),
    params:
        repeats=10,
    script:
        "../scripts/sample_taxons.py"
