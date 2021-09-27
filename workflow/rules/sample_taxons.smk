rule sample_taxons:
    input:
        "resources/speclist.txt"
    log:
        "results/sample_taxons/log.txt"
    output:
        "results/sample_taxons/sample_taxons.txt"
    script:
        "../scripts/sample_taxons.py"
