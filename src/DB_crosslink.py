#######################################################################################
#                               LIPID MAPS REST API                                   # 
#######################################################################################
# El siguiente código:
    # 1. Lee el Excel con una columna ID.
    # 2. Filtra los IDs.
    # 3. Usa la API REST (/rest/compound/lm_id/{lm_id}/all) de LIPID MAPS® 
    #    (https://www.lipidmaps.org/resources/rest) para obtener HMDB, ChEBI, 
    #    InChIKey, etc.
    # 4. Rellena esas columnas (PubChem, KEGG, HMDB, ChEBI, UniProt, InChI, InChIKey) 
    #    en el mismo Excel.
    # 5. Ignora los IDs que no sean de LIPID MAPS, que empiezan por LM.
#######################################################################################

import pandas as pd
import requests
import time
from openpyxl import load_workbook
import re

TARGET_COLS = ["KEGG", "PubChem", "HMDB", "ChEBI", "RefSeq_Id", "UniProt", "InChIKey", "InChI"]

# Definir UniProt regex
uniprot_pattern = r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
    # https://www.uniprot.org/help/accession_numbers

def fetch_lipidmaps_info(lm_id):
    """Fetch cross-references and InChI info from LipidMaps REST API."""
    if not lm_id:
        return None

    # Filtrar IDs
    if lm_id.startswith("LMP"):
        ext = "protein/lmp_id"
    elif re.match(uniprot_pattern, lm_id):
        ext = "protein/uniprot_id"
    else:
        ext = "compound/lm_id"
    url = f"https://www.lipidmaps.org/rest/{ext}/{lm_id}/all"

    # Consultar la base de datos
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        r = requests.get(url, headers=headers, timeout=10)
        r.raise_for_status()
        data = r.json()

        if isinstance(data, list) and not data:
            return None

        # PROTEIN: múltiples filas, usar solo Row1
        if lm_id.startswith("LMP") and isinstance(data, dict) and any(k.startswith('Row') for k in data):
            first_row = data.get('Row1')
            if not first_row:
                return None
            return {
                "LM_ID": lm_id,
                "RefSeq_Id": first_row.get("refseq_id"),
                "UniProt": first_row.get("uniprot_id")
            }    
        # METABOLITE y PROTEIN: dict plano
        elif isinstance(data, dict):
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
        else:
            print(f"[ERROR] Unexpected data format for {lm_id}: {type(data)}")
            return None

    except (requests.exceptions.RequestException, ValueError) as e: 
        print(f"[ERROR] No data for {lm_id} ({e})") 
        return None

# Leer Excel
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
    df.columns = df.iloc[0]     # Primera fila = header
    df = df[2:]                 # Saltamos header y nombre del pathway

    if "ID" not in df.columns:
        print(f"[WARN] No 'ID' column in {sheet_name}, skipped.")
        continue

    # Filtrado y conteo de IDs 
    lm_ids, other_ids, empty = [], [], 0
    for id_raw in df["ID"]:
        id_str = str(id_raw).strip() if pd.notna(id_raw) else '' # si el ID dado no es nulo
        if not id_str:
            empty += 1
        elif id_str.startswith("LM"):
            lm_ids.append(id_str)
        else:
            other_ids.append(id_str)
            
    print(f"Found {len(lm_ids)} LIPID MAPS IDs in {sheet_name}.")
    # print(f"Found {len(other_ids)} IDs wich are not from LIPID MAPS in {sheet_name}:{other_ids}")
    # print(f"Found {empty} IDs wich are empty in {sheet_name}.")

    # Consultar LipidMaps
    results = {}
    for query_id in (lm_ids + other_ids):
        if query_id in results:
            continue
        info = fetch_lipidmaps_info(query_id)
        if info:
            # Guardar mapeo por LM_ID si existe
            if info.get("LM_ID"):
                results[info["LM_ID"]] = info
            # Guardar mapeo por UniProt si existe
            if info.get("UniProt"):
                results[info["UniProt"]] = info       
        time.sleep(0.25)

    # Mapear columnas a índices en la hoja Excel
    header = [cell.value for cell in ws[1]]
    id_idx = header.index("ID") + 1
    target_col_idx = {col: header.index(col) + 1 for col in TARGET_COLS if col in header}

    # Actualizar celdas
    for row in ws.iter_rows(min_row=2, values_only=False):
        id_cell = row[id_idx - 1]
        id_val = str(id_cell.value).strip() if pd.notna(id_cell.value) else ''
        if not id_val:
            continue

        info = results.get(id_val)
        if not info and re.match(uniprot_pattern, id_val):
            info = results.get(id_val)
        
        if info:
            # Reemplazar UniProt → LM_ID si aplica
            if re.match(uniprot_pattern, id_val) and info.get("LM_ID"):
                id_cell.value = info["LM_ID"]
            # Rellenar columnas
            for col_name, col_idx in target_col_idx.items():
                val = info.get(col_name)
                if val:
                    ws.cell(row=id_cell.row, column=col_idx, value=val)

    unmapped_ids = [x for x in other_ids if x not in results]
    print(f"IDs not mapped in {sheet_name}: {unmapped_ids}")

    fin = time.perf_counter()
    total_time.append(fin - inicio)
    print(f"Tiempo de ejecución: {fin - inicio:.2f} s")

# Guardar cambios
output_file = input_file.replace(".xlsx", "_updated.xlsx")
wb.save(output_file)
print(f"Archivo guardado como: {output_file}")
print(f"Tiempo total: {sum(total_time):.2f} s")




#######################################################################################
#                                UniProt REST API                                     
#######################################################################################
# El siguiente código:
    # 1. Lee el Excel con una columna ID.
    # 2. Filtra solo los IDs válidos.
    # 3. Usa la API REST de UniProt (https://www.uniprot.org/help/api_queries) 
    #    para obtener HMDB, ChEBI, 
    #    InChIKey, etc.
    # 4. Rellena esas columnas (PubChem, KEGG, HMDB, ChEBI, UniProt, InChI, InChIKey) 
    #    en el mismo Excel.
    # 5. Ignora los IDs que no sean de LIPID MAPS, que empiezan por LM.



# from UniProtMapper import ProtMapper

# mapper = ProtMapper()

# fields = ["xref_ensembl"]
# result, failed = mapper.get(
#     ids=["Q9NYL5","Q9H2F3","Q9Y6A2"], 
#     fields=fields,
#     )
# print(result)


