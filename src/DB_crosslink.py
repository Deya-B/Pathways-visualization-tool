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


########################################### MAIN ############################################
# Estuctura:
# main.py
# ├── read_input_excel()        # Carga datos desde Excel y devuelve {sheet_name: DataFrame}
# ├── classify_ids()            # Separa LM, UniProt, y otros IDs
# ├── fetch_lipidmaps_info()    # Busca información y cross-referencias de LIPID MAPS
# ├── fetch_pubchem_info()      # Busca información y cross-referencias de PubChem
# ├── fetch_refseq_from_uniprot() # Busca información y cross-referencias de UniProt
# ├── integrate_crossrefs()     # Une resultados de LipidMaps, PubChem y UniProt
# ├── update_excel()            # Escribe los resultados en el Excel
# └── main()                    # Coordina todo

import os
import pandas as pd
import requests
import time
import logging
import re

############################# CONFIGURATION AND CONSTANTS ##############################

# Logging configuration
logging.basicConfig(
    level=logging.INFO, 
    format="%(levelname)s:%(message)s"
    # filename="log.txt"
)

ID_LABEL = "ID"

TARGET_COLS = [         # Columns of interest
    "KEGG", "PubChem", "HMDB", "ChEBI",
    "RefSeq_Id", "UniProt", "InChIKey", "InChI"
]

# UniProt ID Pattern (verified according to official rules)
UNIPROT_PATTERN = r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
    # Ref.: https://www.uniprot.org/help/accession_numbers


############################## CLASSIFICATION FUNCTIONS ################################

def classify_ids(df):
    """Separates LM_IDs, UniProt_IDs and other IDs from the 'ID' column."""
    seen, lm_ids, uniprot_ids, other_ids = set(), [], [], []
    for id_raw in df[ID_LABEL]:
        id_str = str(id_raw).strip() if pd.notna(id_raw) else ""
        if not id_str:
            continue
        for sub in [s.strip() for s in id_str.split(";") if s.strip()]:
            if sub in seen:
                continue
            seen.add(sub)
            if sub.startswith("LM"):
                lm_ids.append(sub)
            elif re.match(UNIPROT_PATTERN, sub):
                uniprot_ids.append(sub)
            else:
                other_ids.append(sub)
    return lm_ids, uniprot_ids, other_ids


############################# FETCHERS (APIs) ##########################################

def fetch_lipidmaps_info(query_id):
    """Search for information using LM_ID or UniProt in the LipidMaps API."""
    if not query_id:
        return None
    if query_id.startswith("LMP"):
        ext = "protein/lmp_id" # extension for LIPID MAPS protein ID queries
    elif re.match(UNIPROT_PATTERN, query_id):
        ext = "protein/uniprot_id" # extension for UniProt protein ID queries
    else:
        ext = "compound/lm_id" # extension for LIPID MAPS metabolite ID queries
    url = f"https://www.lipidmaps.org/rest/{ext}/{query_id}/all"

    # Query the database
    try:
        r = requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=10)
        r.raise_for_status()
        data = r.json()
        if isinstance(data, list) and not data:
            return None
        
        # PROTEIN case: multiple rows, use only Row1
        if isinstance(data, dict) and any(k.startswith('Row') for k in data):
            first_row = data.get('Row1') or {}
            return {
                "LM_ID": first_row.get("lm_id") or first_row.get("lmp_id") or query_id,
                "RefSeq_Id": first_row.get("refseq_id"),
                "UniProt": first_row.get("uniprot_id")
            }
        # METABOLITE/PROTEIN case: flat dict
        elif isinstance(data, dict):
            lm_from_data = data.get("lm_id") or data.get("lmp_id") or None # extract LM id
            return {
                "LM_ID": lm_from_data if lm_from_data else (query_id if query_id.startswith("LM") else None),
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
            logging.error(f"Formato inesperado para {query_id}: {type(data)}")
            return None
    except (requests.exceptions.RequestException, ValueError) as e:
        logging.warning(f"No data for {query_id} ({e})")
        return None


def fetch_pubchem_info(pubchem_id):
    """Busca info por CID en la API de PubChem."""
    try:
        cid = re.sub(r"[^0-9]", "", pubchem_id) # mantener solo dígitos
        if not cid:
            return None
        
        # Extraer InChI e InChIKey
        url_inchi = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChI,InChIKey/JSON"
        r_inchi = requests.get(url_inchi, timeout=10)
        r_inchi.raise_for_status()
        data_inchi = r_inchi.json().get("PropertyTable", {}).get("Properties", [{}])[0]
        # Extraer crossreferences (xrefs)
        url_xrefs = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RegistryID/JSON"
        r_xrefs = requests.get(url_xrefs, timeout=10)
        r_xrefs.raise_for_status()
        data_xrefs = r_xrefs.json().get("InformationList", {}).get("Information", [{}])[0].get("RegistryID", [])
        kegg_id = next((x for x in data_xrefs if re.fullmatch(r"C\d{5}", x)), None)
        chebi_id = next((x for x in data_xrefs if x.startswith("CHEBI:")), None)
        hmdb_id = next((x for x in data_xrefs if x.startswith("HMDB")), None)
        return {
            "PubChem": cid,
            "KEGG": kegg_id,
            "HMDB": hmdb_id,
            "ChEBI": chebi_id,
            "InChI": data_inchi.get("InChI"),
            "InChIKey": data_inchi.get("InChIKey")
        }
    except Exception as e:
        logging.warning(f"PubChem no data for {pubchem_id}: {e}")
        return None


def fetch_refseq_from_uniprot(uniprot_id):
    """Busca acceso RefSeq cruzado en la API de UniProt."""
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}?format=json"
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.json()

        # Extraer IDs RefSeq
        refseq_ids = []
        for xref in data.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "RefSeq":
                refseq_ids.append(xref.get("id"))
        return {
            "UniProt": uniprot_id,
            "RefSeq_Id": ",".join(refseq_ids) if refseq_ids else "NaN"
        }
    except Exception as e:
        logging.warning(f"UniProt no data for {uniprot_id}: {e}")
        return None
    

###################### INTEGRACIÓN DE CROSS-REFERENCIAS #########################

def integrate_crossrefs(df):
    """Une resultados de LipidMaps, PubChem y UniProt para todos los IDs."""
    # Clasificar IDs
    seen_queries = set()
    lm_ids, uniprot_ids, other_ids = classify_ids(df)
    # Totales por tipo:
    logging.info(f"  LIPID MAPS: {len(lm_ids)}  |  UniProt: {len(uniprot_ids)}  |  Otros: {len(other_ids)}")

    results = {}

    # Consultar LipidMaps
    for query_id in lm_ids + other_ids:
        info = fetch_lipidmaps_info(query_id)
        seen_queries.add(query_id)
        if not info:
            # time.sleep(0.25)
            continue
        lm_key = info.get("LM_ID")
        if lm_key:
            results[lm_key] = info
        uni_field = info.get("UniProt")
        if uni_field:
            for uni in [u.strip() for u in str(uni_field).split(";") if u.strip()]:
                results[uni] = info
        results[query_id] = info
        # time.sleep(0.25)
    # IDs no crossreferenciadas
    not_crossreferenced = set()
    for id_raw in df[ID_LABEL]:
        id_str = str(id_raw).strip() if pd.notna(id_raw) else ''
        if not id_str:
            continue
        for sub in [s.strip() for s in id_str.split(";") if s.strip()]:
            if sub not in results:
                not_crossreferenced.add(sub)

    # Consultar UniProt para ids sin cruzar
    for uid in uniprot_ids:
        if uid not in results:
            info = fetch_refseq_from_uniprot(uid)
            seen_queries.add(uid)
            if info:
                results[uid] = info
            else:
                not_crossreferenced.add(uid)
            # time.sleep(0.2)

    # Consultar PubChem para ids numéricos sin cruzar
    for qid in sorted(not_crossreferenced):
        if re.match(r"^\d+$", qid):
            info = fetch_pubchem_info(qid)
            seen_queries.add(qid)
            if info:
                results[qid] = info
                not_crossreferenced.discard(qid)
            # time.sleep(0.3)

    # Excluir UniProt del resumen
    not_crossreferenced_final = [
        x for x in sorted(seen_queries)
        if x not in results and not re.match(UNIPROT_PATTERN, x)
    ]
    
    return results, not_crossreferenced_final
    

############################# ESCRITURA en DataFrame #############################

def update_df_with_crossrefs(df, results):
    """Actualiza el DataFrame con los crossrefs encontrados (modifica las TARGET_COLS)."""
    for idx, row in df.iterrows():
        id_str = str(row[ID_LABEL]).strip() if pd.notna(row[ID_LABEL]) else ""
        sub_ids = [x.strip() for x in id_str.split(";") if x.strip()]
        for col in TARGET_COLS:
            vals = []
            for sub_id in sub_ids:
                info = results.get(sub_id)
                if info and info.get(col):
                    vals.append(str(info.get(col)))
            if vals:
                df.at[idx, col] = ";".join(vals)
            else:
                df.at[idx, col] = ""
    return df


########################## MAIN CONTROLADOR #####################################

def main(input_file):
    # Read input
    total_time = []
    start = time.perf_counter()
    df = pd.read_csv(input_file, sep='\t', encoding="cp1252").dropna(axis=0, how='all').dropna(axis=1, how='all')
                                            # con ".dropna" quitamos las filas y columnas vacías
    base_filename = os.path.basename(input_file)
    logging.info(f" --- Procesando archivo: {base_filename} ---")
    
    # Busqueda de crossreferencias
    results, not_crossreferenced_final = integrate_crossrefs(df)
    df = update_df_with_crossrefs(df, results)
    
    # Guardar el df resultante en tsv
    output_folder = "c:/Users/dborrotoa/Desktop/TFM/pathways_updated/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    base_filename_upd = base_filename.replace(".txt", "_updated.txt")
    output_path = os.path.join(output_folder,base_filename_upd)
    df.to_csv(output_path, sep='\t')

    # Tiempo de ejecución y registro de la ejecución
    end = time.perf_counter()
    total_time.append(end - start)
    logging.info(f"     Tiempo ejecución hoja: {end - start:.2f} s")
    logging.info(f"     IDs no crossreferenciados: {len(not_crossreferenced_final)}")
    logging.info(f" Archivo guardado en: {output_path}\n")


############################ ENTRY POINT ########################################

if __name__ == "__main__":
    pathways_folder = "c:/Users/dborrotoa/Desktop/TFM/pathways_raw/"
    tsv_files = [f for f in os.listdir(pathways_folder) if f.endswith('.txt') or f.endswith('.tsv')]
    for file in tsv_files:
        INPUT_FILE = os.path.join(pathways_folder, file)
        main(INPUT_FILE)

