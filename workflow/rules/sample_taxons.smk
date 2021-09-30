rule sample_taxons:
    input:
        "resources/speclist.txt"
    log:
        "logs/sample_taxons.txt"
    output:
        "results/sample_taxons/sample_taxons.txt"
    script:
        "../scripts/sample_taxons.py"
