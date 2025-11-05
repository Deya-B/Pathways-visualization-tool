# Fetch other databases (HMDB, ChEBI, Ensembl, UniProt) and the InChI_key through the LM_ID
# by using LIPID MAPS® REST service (https://www.lipidmaps.org/resources/rest)

# PARTS:
# 1. invariant base URL
base_URL = "https://www.lipidmaps.org/rest"

# 2. INPUT 
# This specification is composed of 3 required parameters separated by forward slashes.
# The first parameter is the context, either "compound", "gene" or "protein", 
# each of which has a separate list of input items associated with it (2nd parameter). 
# The 3rd parameter is an appropriate input value for the chosen item.

# 3. OUTPUT
# Described in url




import requests

def fetch_lipidmaps_info(lm_id):
    url = f"https://www.lipidmaps.org/rest/lmwebservice.php?LM_ID={lm_id}"
    r = requests.get(url)
    r.raise_for_status()
    data = r.json()
    return data

lm_id = "LMFA01010001"  # ejemplo: palmitic acid
info = fetch_lipidmaps_info(lm_id)

print(info.keys())
# Este devuelve un JSON con todos los campos disponibles.



# Para EXTRAER campos específicos
def extract_crossrefs(data):
    refs = data.get("EXTERNAL_REFERENCES", {})
    return {
        "HMDB": refs.get("HMDB_ID"),
        "ChEBI": refs.get("CHEBI_ID"),
        "Ensembl": refs.get("ENSEMBL_ID"),
        "UniProt": refs.get("UNIPROT_ID"),
        "InChIKey": data.get("INCHI_KEY")
    }

refs = extract_crossrefs(info)
print(refs)


# Integrar varios LM_IDs
import pandas as pd

lm_ids = ["LMFA01010001", "LMFA02000002", "LMFA03000003"]
results = []

for lm_id in lm_ids:
    data = fetch_lipidmaps_info(lm_id)
    refs = extract_crossrefs(data)
    refs["LM_ID"] = lm_id
    results.append(refs)

df = pd.DataFrame(results)
print(df)



# Para guardar los resultados en CSV 
df.to_csv("lipidmaps_crossrefs.csv", index=False)