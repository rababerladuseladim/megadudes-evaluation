from snakemake.remote.HTTP import RemoteProvider
HTTP = RemoteProvider()

rule sample_taxons:
    input:
        HTTP.remote("https://www.uniprot.org/docs/speclist.txt")
    output:
        "results/sample_taxons/sample_taxons.txt"
    script:
        "../scripts/sample_taxons.py"
