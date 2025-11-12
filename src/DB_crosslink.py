
#######################################################################################
#                 LIPID MAPS REST API + PubChem REST API                              # 
#######################################################################################
# El siguiente código:
    # 1. Lee el Excel con los ID's de LIPID MAPS y PubChem.
    # 2. Filtra los IDs: de LIPID MAPS para metabolitos o proteínas, y los de PubChem. 
    # 3. 
    #  a) Usa la API REST de LIPID MAPS® (/rest/compound/lm_id/{lm_id}/all)  
    #       (https://www.lipidmaps.org/resources/rest) para obtener la InChI e InChIKey de
    #       cada metabolito y las cross-refs de los IDs con KEGG, PubChem, HMDB y ChEBI
    #       (metabolitos) y RefSeq_ID y UniProt en el caso de proteínas.
    #  b) Usa las API REST de PubChem (https://pubchem.ncbi.nlm.nih.gov/docs/rdf-rest) 
    #       para obtener, tanto los datos de la InChI e InChIKey de cada metabolito, 
    #       como las cross-refs de los IDs con KEGG, HMDB y ChEBI
    # 4. Genera una copia del Excel actualizado "_updated.xlsx".
#######################################################################################

import pandas as pd
import requests
import time
from openpyxl import load_workbook
import re

def fetch_pubchem_info(pubchem_id):
    """Fetch compound info from PubChem REST API (using CID)."""
    try:
        cid = re.sub(r"[^0-9]", "", pubchem_id)  # mantener solo dígitos
        if not cid:
            return None
        
        # Extraer InChI e InChIKey
        url_inchi = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChI,InChIKey/JSON"
        r_inchi = requests.get(url_inchi, timeout=10)
        r_inchi.raise_for_status()
        data_inchi = (
            r_inchi.json()
            .get("PropertyTable", {})
            .get("Properties", [{}])[0]
        )
        # Extraer crossreferences (xrefs)
        url_xrefs = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RegistryID/JSON"
        r_xrefs = requests.get(url_xrefs, timeout=10)
        r_xrefs.raise_for_status()
        data_xrefs = (
            r_xrefs.json()
            .get("InformationList", {})
            .get("Information", [{}])[0]
            .get("RegistryID", [])
        )

        kegg_id = next((x for x in data_xrefs if re.fullmatch(r"C\d{5}", x)), None) # C+numeros
        chebi_id = next((x for x in data_xrefs if x.startswith("CHEBI:")), None)
        hmdb_id = next((x for x in data_xrefs if x.startswith("HMDB")), None)

        res = {
            "PubChem": cid,
            "KEGG": kegg_id,
            "HMDB": hmdb_id,
            "ChEBI": chebi_id,
            "InChI": data_inchi.get("InChI"),
            "InChIKey": data_inchi.get("InChIKey")
        }
        return res
    except Exception as e:
        print(f"[WARN] PubChem no data for {pubchem_id}: {e}")
        return None
###### NEW

# TODO: Allow for capital/lowercase letters in column names
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
            print(f"[ERROR] Unexpected data format for {query_id}: {type(data)}")
            return None

    except (requests.exceptions.RequestException, ValueError) as e: 
        print(f"[ERROR] No data for {query_id} ({e})") 
        return None

# Leer Excel
input_file = "c:/Users/dborrotoa/Desktop/TFM/PathwayBA_list.xlsx"
# input_file = "C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/PathwayBA_list.xlsx"
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
    # Intentar recuperar información de PubChem
    external_results = {}
    for unknown_id in sorted(not_crossreferenced):
        if re.match(r"^\d+$", unknown_id):  # numérico → PubChem CID
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

output_file = input_file.replace(".xlsx", "_updated.xlsx")
wb.save(output_file)
print(f"Archivo guardado como: {output_file}")
print(f"Tiempo total: {sum(total_time):.2f} s")
