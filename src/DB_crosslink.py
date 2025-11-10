# El siguiente código:
    # 1. Lee el Excel con una columna ID.
    # 2. Filtra solo los IDs válidos de LIPID MAPS (LM...).
    # 3. Usa la API REST (/rest/compound/lm_id/{lm_id}/all) de LIPID MAPS® (https://www.lipidmaps.org/resources/rest)
    #    para obtener HMDB, ChEBI, InChIKey, etc.
    # 4. Rellena esas columnas (PubChem, KEGG, HMDB, ChEBI, UniProt, InChI, InChIKey) en el mismo Excel.
    # 5. Ignora los IDs que no sean de LIPID MAPS, que empiezan por LM.
####################################################################################################################

import pandas as pd
import requests
import time
from openpyxl import load_workbook

TARGET_COLS = ["KEGG", "PubChem", "HMDB", "ChEBI", "RefSeq_Id", "UniProt", "InChIKey", "InChI"]

def fetch_lipidmaps_info(lm_id):
    """Fetch cross-references and InChI info from LipidMaps REST API."""

    proteins_ext = "protein/lmp_id"
    metabolites_ext = "compound/lm_id"
    if lm_id.startswith("LMP"):
        ext = proteins_ext
    else:
        ext = metabolites_ext

    url = f"https://www.lipidmaps.org/rest/{ext}/{lm_id}/all"

    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        r = requests.get(url, headers=headers, timeout=10)
        r.raise_for_status()
        data = r.json()

        if isinstance(data, list) and not data:
            return None

        if not isinstance(data, dict):
            print(f"[ERROR] Unexpected data format for {lm_id}: {type(data)}")
            return None

        return {
            "LM_ID": lm_id,
            "KEGG": data.get("kegg_id"),
            "PubChem": data.get("pubchem_cid"),
            "HMDB": data.get("hmdb_id"),
            "ChEBI": data.get("chebi_id"),
            "RefSeq_Id": data.get("refseq_id"),
            "UniProt": data.get("uniprot_id"),
            "InChIKey": data.get("inchi_key"),
            "InChI": data.get("inchi")
        }
    except (requests.exceptions.RequestException, ValueError): 
        print(f"[ERROR] No data for {lm_id}") 
        return None

# 1. Leer Excel
# input_file = "c:/Users/dborrotoa/Desktop/TFM/PathwayBA_list.xlsx"
input_file = "C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/PathwayBA_list.xlsx"
wb = load_workbook(input_file)
total_time = []

for sheet_name in wb.sheetnames:
    inicio = time.perf_counter()
    ws = wb[sheet_name]
    print(f"\n--- Processing sheet: {sheet_name} ---")

    # Cargar datos de la hoja con pandas para trabajar más cómodo
    df = pd.DataFrame(ws.values)
    df.columns = df.iloc[0] # Primera fila = header
    df = df[2:] # Saltamos header y nombre del pathway

    if "ID" not in df.columns:
        print(f"[WARN] No 'ID' column in {sheet_name}, skipped.")
        continue

    # Filtrado y conteo de IDs 
    lm_ids = []
    other_ids = []
    empty = 0
    for id_raw in df["ID"]:
        id_str = str(id_raw).strip() if pd.notna(id_raw) else '' # si el ID dado no es nulo
        if not id_str:
            empty += 1
        elif id_str.startswith("LM"):
            lm_ids.append(id_str)
        else:
            other_ids.append(id_str)
            
    print(f"Found {len(lm_ids)} LIPID MAPS IDs in {sheet_name}.")
    print(f"Found {len(other_ids)} IDs wich are not from LIPID MAPS in {sheet_name}.")
    print(f"Found {empty} IDs wich are empty in {sheet_name}.")

    # Consultar LipidMaps
    results = {}
    for lm_id in lm_ids:
        info = fetch_lipidmaps_info(lm_id)
        if info:
            results[lm_id] = info
        time.sleep(0.3)

    # Mapear columnas a índices en la hoja Excel
    header = [cell.value for cell in ws[1]]
    col_idx = {col: header.index(col) + 1 for col in TARGET_COLS if col in header}
    id_idx = header.index("ID") + 1

    # Actualizar celdas directamente
    for row in ws.iter_rows(min_row=2, values_only=False):
        lm_id = str(row[id_idx - 1].value)
        if lm_id in results:
            for col_name, value in results[lm_id].items():
                if col_name in col_idx and value is not None:
                    ws.cell(row=row[0].row, column=col_idx[col_name], value=value)

    fin = time.perf_counter()
    total_time.append(fin - inicio)
    print(f"Tiempo de ejecución: {fin - inicio:.2f} s")

# Guardar cambios
output_file = input_file.replace(".xlsx", "_updated.xlsx")
wb.save(output_file)
print(f"Archivo guardado como: {output_file}")
print(f"Tiempo total: {sum(total_time):.2f} s")





# import requests, json

# url = "https://www.lipidmaps.org/rest/protein/lmp_id/LMP002143/all"
# print(json.dumps(requests.get(url).json(), indent=2))