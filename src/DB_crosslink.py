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

def fetch_lipidmaps_info(query_id):
    """Fetch cross-references and InChI info from LipidMaps REST API."""
    if not query_id:
        return None

    # Filtrar IDs
    if query_id.startswith("LMP"):
        ext = "protein/lmp_id"
    elif re.match(uniprot_pattern, query_id):
        ext = "protein/uniprot_id"
    else:
        ext = "compound/lm_id"
    url = f"https://www.lipidmaps.org/rest/{ext}/{query_id}/all"

    # Consultar la base de datos
    try:
        r = requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=10)
        r.raise_for_status()
        data = r.json()
        if isinstance(data, list) and not data:
            return None

        # PROTEIN: múltiples filas, usar solo Row1
        if isinstance(data, dict) and any(k.startswith('Row') for k in data):
            first_row = data.get('Row1') or {}
            return {
                "LM_ID": first_row.get("lm_id") or first_row.get("lmp_id") or query_id,
                "RefSeq_Id": first_row.get("refseq_id"),
                "UniProt": first_row.get("uniprot_id")
            }    
        # METABOLITE y PROTEIN: dict plano
        elif isinstance(data, dict):
            # extraer LM id
            lm_from_data = data.get("lm_id") or data.get("lmp_id") or None
            return {
                "LM_ID": lm_from_data if lm_from_data 
                            else (query_id if query_id.startswith("LM") else None),
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

    # Filtrado y conteo de IDs 
    seen_queries = set()
    lm_ids, other_ids = [], []
    for id_raw in df["ID"]:
        id_str = str(id_raw).strip() if pd.notna(id_raw) else '' # si el ID dado no es nulo
        if not id_str:
            continue
        for sub in [s.strip() for s in id_str.split(";") if s.strip()]:
            if sub in seen_queries:
                continue
            seen_queries.add(sub)
            if sub.startswith("LM"):
                lm_ids.append(sub)
            else:
                other_ids.append(sub)


    print(f"Found {len(lm_ids)} LIPID MAPS IDs in {sheet_name}.")
    print(f"Found {len(other_ids)} IDs wich are not from LIPID MAPS in {sheet_name}:{other_ids}")

    # Consultar LipidMaps  y crear index de búsqueda
    results = {}
    for query_id in (lm_ids + other_ids):
        info = fetch_lipidmaps_info(query_id)
        if not info:
            time.sleep(0.25)
            continue

        # LM ID
        lm_key = info.get("LM_ID")
        if lm_key:
            results[lm_key] = info

        # UniProt ID
        uni_field = info.get("UniProt")
        if uni_field:
            for uni in [u.strip() for u in str(uni_field).split(";") if u.strip()]:
                results[uni] = info
        
        # Indexar por la consulta original (por si no hay LM o UniProt)
        results[query_id] = info
        time.sleep(0.25)

    # Cálculo de IDs no crossreferenciados
    not_crossreferenced = set()
    for id_raw in df["ID"]:
        id_str = str(id_raw).strip() if pd.notna(id_raw) else ''
        if not id_str:
            continue
        for sub in [s.strip() for s in id_str.split(";") if s.strip()]:
            if sub not in results:
                not_crossreferenced.add(sub)
    
    if not_crossreferenced:
        print(f"--- [{len(not_crossreferenced)}] IDs finished with no crossreference in {sheet_name}: {sorted(not_crossreferenced)}")
    else:
        print(f"Todas las IDs cruzadas en {sheet_name}.")


###### NEW
    # Intentar recuperar información de PubChem y ChEBI
    external_results = {}
    for unknown_id in sorted(not_crossreferenced):
        if unknown_id.startswith("CHEBI:"):
            info = fetch_chebi_info(unknown_id)
        elif re.match(r"^\d+$", unknown_id):  # numérico → PubChem CID
            info = fetch_pubchem_info(unknown_id)
        else:
            continue

        if info:
            external_results[unknown_id] = info
        time.sleep(0.3)

    # Unir resultados externos al diccionario principal
    results.update(external_results)
###### NEW



    # Mapear columnas a índices en la hoja Excel
    header = [cell.value for cell in ws[1]]
    if "ID" not in header:
        print(f"[WARN] No 'ID' header in sheet {sheet_name}, skip.")
        continue
    id_idx = header.index("ID") + 1
    target_col_idx = {col: header.index(col) + 1 for col in TARGET_COLS if col in header}

    # Actualizar celdas
    for row in ws.iter_rows(min_row=2, values_only=False):
        id_cell = row[id_idx - 1]
        id_val = str(id_cell.value).strip() if pd.notna(id_cell.value) else ''
        if not id_val:
            continue

        sub_ids = [x.strip() for x in id_val.split(";") if x.strip()]
        merged_info = {k: [] for k in TARGET_COLS}
        new_id_val = []

        for sub_id in sub_ids:
            info = results.get(sub_id)
            if not info and re.match(uniprot_pattern, sub_id):
                info = results.get(sub_id)
            if info:
                # reemplazar UniProt por LM si existe LM_ID real
                lm_real = info.get("LM_ID")
                if lm_real:
                    new_id_val.append(lm_real)
                else:
                    new_id_val.append(sub_id)

                # para cada campo target, añadir valores evitando duplicados
                for col in TARGET_COLS:
                    v = info.get(col)
                    if v:
                        # si la API devuelve varios dentro del mismo campo, separar por ';'
                        for part in [p.strip() for p in str(v).split(";") if p.strip()]:
                            if part not in merged_info[col]:
                                merged_info[col].append(part)
            else:
                # si no hay info, conservar el sub-id sin cambios
                new_id_val.append(sub_id)

        # escribir nuevo ID concatenado (mantener orden encontrado)
        id_cell.value = ";".join(new_id_val)

        # escribir columnas objetivo concatenadas
        for col_name, col_idx in target_col_idx.items():
            vals = merged_info.get(col_name) or []
            if vals:
                ws.cell(row=id_cell.row, column=col_idx, value=";".join(vals))

    fin = time.perf_counter()
    total_time.append(fin - inicio)
    print(f"Tiempo de ejecución: {fin - inicio:.2f} s")

output_file = input_file.replace(".xlsx", "_updated2.xlsx")
wb.save(output_file)
print(f"Archivo guardado como: {output_file}")
print(f"Tiempo total: {sum(total_time):.2f} s")


#######################################################################################
#                           CheBI + PubChem REST API                                     
#######################################################################################
# El siguiente código:
    # 1. Con los IDs no crossreferenciados que pertenezcan a CheBi o PubChem del
    #    not_crossreferenced = set()
    # 2. Usa las API REST (https://pubchem.ncbi.nlm.nih.gov/docs/rdf-rest) 
    #   (https://jcheminf.biomedcentral.com/articles/10.1186/s13321-016-0123-9)
    #    para obtener los datos
    # 4. Rellena las columnas (PubChem, KEGG, HMDB, ChEBI, UniProt, InChI, InChIKey) 
    #    del Excel

# --- Processing sheet: Sheet1 ---
# Todas las IDs cruzadas en Sheet1.

# --- Processing sheet: Sheet2 ---
# IDs not crossreferenced in Sheet2: 
# ['155920193', '42776911', 'A7AZH2', 'A7B3K3', 'A7B4V1', 
# 'CHEBI:133678', 'CHEBI:172393', 'P0AET8', 'P21215']

# --- Processing sheet: Sheet3 ---
# IDs not crossreferenced in Sheet3: 
# ['71448933', '91828270', 'A7AZH2', 'A7B3K3', 'A7B4V1', 'C8WGQ3', 
# 'C8WMP0', 'P0AET8', 'Q5LA59']

# --- Processing sheet: Sheet4 ---
# IDs not crossreferenced in Sheet4: 
# ['137333454', 'A7B4V1', 'P0AET8', 'P52843', 'Q4LDG0', 'Q91W64']

# --- Processing sheet: Sheet5 ---
# IDs not crossreferenced in Sheet5: 
# ['11966205', '134364', '644071', 'A7AZH2', 'A7B3K3', 'B0NAQ4', 'C8WGQ3', 
# 'C8WMP0', 'P07914', 'P19337', 'P19409', 'P19410', 'P19412', 'P19413', 'P32370']

# https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/644071/property/InChI,InChIKey,Title/JSON
# https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/644071/xrefs/RegistryID/JSON


###### NEW
def fetch_chebi_info(chebi_id):
    """Fetch compound info from ChEBI REST API."""
    try:
        chebi_id = chebi_id.replace(" ", "")
        if not chebi_id.startswith("CHEBI:"):
            return None
        url = f"https://www.ebi.ac.uk/chebi/ws/rest/{chebi_id}?format=json"
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.json()
        res = {
            "ChEBI": chebi_id,
            "KEGG": None,
            "PubChem": None,
            "HMDB": None,
            "InChI": data.get("inchi"),
            "InChIKey": data.get("inchiKey")
        }
        return res
    except Exception as e:
        print(f"[WARN] ChEBI no data for {chebi_id}: {e}")
        return None


def fetch_pubchem_info(pubchem_id):
    """Fetch compound info from PubChem REST API (using CID)."""
    try:
        cid = re.sub(r"[^0-9]", "", pubchem_id)  # mantener solo dígitos
        if not cid:
            return None
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChI,InChIKey,Title/JSON"
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.json().get("PropertyTable", {}).get("Properties", [{}])[0]
        res = {
            "PubChem": cid,
            "KEGG": None,
            "HMDB": None,
            "ChEBI": None,
            "InChI": data.get("InChI"),
            "InChIKey": data.get("InChIKey")
        }
        return res
    except Exception as e:
        print(f"[WARN] PubChem no data for {pubchem_id}: {e}")
        return None
###### NEW