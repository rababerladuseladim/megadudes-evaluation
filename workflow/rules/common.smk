from snakemake.utils import validate
import pandas as pd


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


##### load config and sample sheets #####


configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_table("config/samples.tsv", dtype=str)
    .set_index("sample_name", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")
