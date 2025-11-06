# El siguiente código:
    # 1. Lee el Excel con una columna ID.
    # 2. Filtra solo los IDs válidos de LIPID MAPS (LM...).
    # 3. Usa la API REST (/rest/compound/lm_id/{lm_id}/all) de LIPID MAPS® (https://www.lipidmaps.org/resources/rest)
    #    para obtener HMDB, ChEBI, InChIKey, etc.
    # 4. Rellena esas columnas (HMDB, ChEBI, Ensembl, UniProt, InChI) en el mismo Excel.
    # 5. Ignora los IDs que no sean de LIPID MAPS, que empiezan por LM.
####################################################################################################################

import pandas as pd
import requests
from time import sleep

def fetch_lipidmaps_info(lm_id):
    """Consulta LIPID MAPS REST API para un LM_ID y devuelve los campos relevantes."""
    url = f"https://www.lipidmaps.org/rest/compound/lm_id/{lm_id}/all"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.json()
        return {
            "LM_ID": lm_id,
            "Name": data.get("name"),
            "HMDB": data.get("hmdb_id"),
            "ChEBI": data.get("chebi_id"),
            "KEGG": data.get("kegg_id"),
            "PubChem": data.get("pubchem_cid"),
            "InChI": data.get("inchi"),
            "InChIKey": data.get("inchi_key")
        }
    except (requests.exceptions.RequestException, ValueError):
        print(f"[ERROR] No data for {lm_id}")
        return None
    

# 1. Leer Excel
input_file = "c:/Users/dborrotoa/Desktop/TFM/PathwayBA_list.xlsx"
output_file = "c:/Users/dborrotoa/Desktop/TFM/PathwayBA_crossrefs.xlsx"

xls = pd.ExcelFile(input_file)
print(f"Sheets found: {xls.sheet_names}")

writer = pd.ExcelWriter(output_file, engine="openpyxl")


# 2. Procesar cada hoja
for sheet_name in xls.sheet_names:
    print(f"\n--- Processing sheet: {sheet_name} ---")
    df = pd.read_excel(xls, sheet_name)

    # verificar que exista columna ID
    if "ID" not in df.columns:
        print(f"[WARN] No 'ID' column in {sheet_name}, skipped.")
        df.to_excel(writer, sheet_name=sheet_name, index=False)
        continue

    # 3. Filtrar IDs válidos (LIPID MAPS)
    lm_ids = [x for x in df["ID"].astype(str) if x.startswith("LM")]
    print(f"Found {len(lm_ids)} LIPID MAPS IDs in {sheet_name}.")

    results = []
    for lm_id in lm_ids:
        info = fetch_lipidmaps_info(lm_id)
        if info:
            results.append(info)
        sleep(0.3)

    if results:
        crossrefs_df = pd.DataFrame(results)
        merged = df.merge(crossrefs_df, how="left", left_on="ID", right_on="LM_ID")

        # combinar columnas duplicadas si existieran
        for col in ["HMDB", "ChEBI", "InChI"]:
            if f"{col}_y" in merged.columns and f"{col}_x" in merged.columns:
                merged[col] = merged[f"{col}_y"].combine_first(merged[f"{col}_x"])

        merged = merged.drop(columns=[c for c in merged.columns if c.endswith(("_x", "_y")) or c == "LM_ID"])
    else:
        merged = df

    # 4. Guardar hoja procesada en el Excel de salida
    merged.to_excel(writer, sheet_name=sheet_name, index=False)

writer.close()
print(f"\nFinished. Results saved to: {output_file}")