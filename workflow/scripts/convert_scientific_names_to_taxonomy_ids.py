import requests
import sys

LOG_HANDLE = sys.stderr


class TaxonomyConnector:
    url = "https://www.ebi.ac.uk/proteins/api/taxonomy/name/"

    def __init__(self):
        self.session = requests.Session()

    def get_taxonomy_id(self, scientific_name: str) -> int | None:
        params = {
            "pageNumber": "1",
            "pageSize": "100",
            "searchType": "EQUALSTO",
            "fieldName": "SCIENTIFICNAME",
        }
        url = self.url + scientific_name.replace(" ", "%20")
        response = self.session.get(url, params=params)
        # ebi api returns 404 if search result is empty
        if response.status_code == 404:
            return None
        response.raise_for_status()
        taxonomies = response.json()["taxonomies"]
        if len(taxonomies) != 1:
            return None
        return taxonomies[0]["taxonomyId"]


def test_get_taxonomy_id():
    taxonomy = TaxonomyConnector()
    assert taxonomy.get_taxonomy_id("Actinomyces sp. oral taxon 448") == 712124
    assert taxonomy.get_taxonomy_id("foo") is None


def convert_species_names_to_taxids(species_file, outfile):
    taxonomy = TaxonomyConnector()
    taxonomy_ids = []
    with open(species_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            if tax_id := taxonomy.get_taxonomy_id(line.replace(" sp ", " sp. ")):
                taxonomy_ids.append(tax_id)
            else:
                print(f"No hit for: line", file=LOG_HANDLE)
    with open(outfile, "w") as fh:
        fh.write("\n".join(map(str, taxonomy_ids)) + "\n")


def test_convert_human_microbiome_species_to_taxonomy_id(tmpdir):
    outfile = tmpdir / "wasd"
    species_file = "resources/human_microbiome_species.txt"
    convert_species_names_to_taxids(species_file=species_file, outfile=outfile)


if snakemake := globals().get("snakemake"):
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        convert_species_names_to_taxids(
            species_file=snakemake.input[0],
            outfile=snakemake.output[0],
        )
