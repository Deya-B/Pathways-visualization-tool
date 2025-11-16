#######################################################################################
#                      LIPID MAPS/PubChem/UniProt REST APIs                           # 
#######################################################################################
# El siguiente código:
    # 1. Lee el Excel con los ID's de LIPID MAPS, PubChem y UniProt.
    # 2. Filtra los IDs. 
    # 3. 
    #  a) Usa la API REST de LIPID MAPS® (https://www.lipidmaps.org/resources/rest) 
    #       desde 2 different endpoints para obtener la InChI e InChIKey de
    #       cada metabolito y las cross-refs de los IDs con KEGG, PubChem, HMDB y ChEBI
    #       (metabolitos) y RefSeq_ID y UniProt en el caso de proteínas.
    #  b) Usa las API REST de PubChem (https://pubchem.ncbi.nlm.nih.gov/docs/rdf-rest) 
    #       para obtener, tanto los datos de la InChI e InChIKey de cada metabolito, 
    #       como las cross-refs de los IDs con KEGG, HMDB y ChEBI
    # 4. Genera una copia del Excel actualizado "_updated.xlsx".
#######################################################################################


########################################### MAIN ######################################
# Estuctura:
# main.py
# ├── classify_ids()            # Separa LM, UniProt, y otros IDs
# ├── fetch_lipidmaps_info()    # Busca información y cross-referencias de LIPID MAPS
# ├── fetch_pubchem_info()      # Busca información y cross-referencias de PubChem
# ├── fetch_refseq_from_uniprot() # Busca información y cross-referencias de UniProt
# ├── integrate_crossrefs()     # Une resultados de LipidMaps, PubChem y UniProt
# ├── update_excel()            # Escribe los resultados en el Excel
# └── main()                    # Coordina todo

import os
import pandas as pd
import numpy as np
import requests
import time
import logging
import re

############################# CONFIGURATION AND CONSTANTS #############################

# Folders
# INPUT_FOLDER = "c:/Users/dborrotoa/Desktop/TFM/pathways_raw/"
INPUT_FOLDER = "C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/pathways_raw"
# OUTPUT_FOLDER = "c:/Users/dborrotoa/Desktop/TFM/pathways_updated/"
OUTPUT_FOLDER = "C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/pathways_updated"

# Logging configuration
logging.basicConfig(
    level=logging.INFO, 
    format="%(levelname)s:%(message)s"
    # filename="log.txt"
)

# Variables with accepted variants
PCHEM = ["pubchem cid", "pubchem"]
LM = ["lipidmaps"]
UPROT = ["uniprot"]
# TODO: Columns
# COL_ORDER = [header]

# UniProt ID Pattern (verified according to official rules)
UNIPROT_PATTERN = r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
    # Ref.: https://www.uniprot.org/help/accession_numbers


############################## CLASSIFICATION FUNCTIONS ###############################

def classify_ids(ids_db_list):
    """Filtering the 'ID' column per DataBase."""
    # Extract ID's (column 0) where DataBase is matching
    lm_ids = ids_db_list[np.isin(ids_db_list[:, 1], LM), 0]
    uni_ids =  ids_db_list[np.isin(ids_db_list[:, 1], UPROT), 0]
    pchem_cids = ids_db_list[np.isin(ids_db_list[:, 1], PCHEM), 0]
    return lm_ids,uni_ids,pchem_cids

############################# FETCHERS (APIs) #########################################

def fetch_lipidmaps_info(query_id):
    """Search for information using LipidMaps/UniProt IDs in the LipidMaps REST API."""
    # Obtain context and identifier type segment for DataBase endpoint
        # https://www.lipidmaps.org/rest/{context}/{input_item}/{input_value}/{output_item}
    if query_id.startswith("LMP"):      # /{context}/{input_item}
        endpoint = "protein/lmp_id"     # for LM protein ID queries
    elif re.match(UNIPROT_PATTERN, query_id):
        endpoint = "protein/uniprot_id" # for UniProt protein ID queries
    else:
        endpoint = "compound/lm_id"     # for LipidMaps metabolite ID queries
    
    url = f"https://www.lipidmaps.org/rest/{endpoint}/{query_id}/all"

    try: # Query the database
        r = requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=10)
        r.raise_for_status()
        data = r.json()
        if isinstance(data, list) and not data:
            return None
        
        # PROTEIN case: multiple rows in response(r), use only Row1
        if isinstance(data, dict) and any(k.startswith('Row') for k in data):
            first_row = data.get('Row1') or {}
            return {
                # "LM_ID": first_row.get("lm_id") or first_row.get("lmp_id") or query_id,
                "RefSeq_Id": first_row.get("refseq_id"),
                "UniProt": first_row.get("uniprot_id")
            }
        # METABOLITE/PROTEIN case: flat dict in response
        elif isinstance(data, dict):
            # lm_from_data = data.get("lm_id") or data.get("lmp_id") or None # extract LM id
            return {
                # "LM_ID": lm_from_data if lm_from_data else (query_id if query_id.startswith("LM") else None),
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
            logging.error(f"Unexpected format for {query_id}: {type(data)}")
            return None
    except (requests.exceptions.RequestException, ValueError) as e:
        logging.warning(f"No data for {query_id} ({e})")
        return None


def fetch_refseq_from_uniprot(query_id):
    """Looking for cross-RefSeq access in the UniProt API."""
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{query_id}?format=json"
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.json()

        # Extraer IDs RefSeq
        refseq_ids = []
        for xref in data.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "RefSeq":
                refseq_ids.append(xref.get("id"))
        return {
            "UniProt": query_id,
            "RefSeq_Id": ",".join(refseq_ids) if refseq_ids else "NaN"
        }
    except Exception as e:
        logging.warning(f"UniProt no data for {query_id}: {e}")
        return None
    

def fetch_pubchem_info(query_id):
    """Search for information by CID in the PubChem API.""" 
    try:    
        # Extract InChI and InChIKey
        url_inchi = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{query_id}/property/InChI,InChIKey/JSON"
        r_inchi = requests.get(url_inchi, timeout=10)
        r_inchi.raise_for_status()
        data_inchi = r_inchi.json().get("PropertyTable", {}).get("Properties", [{}])[0]
        # Extract crossreferences (xrefs)
        url_xrefs = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{query_id}/xrefs/RegistryID/JSON"
        r_xrefs = requests.get(url_xrefs, timeout=10)
        r_xrefs.raise_for_status()
        data_xrefs = r_xrefs.json().get("InformationList", {}).get("Information", [{}])[0].get("RegistryID", [])
        kegg_id = next((x for x in data_xrefs if re.fullmatch(r"C\d{5}", x)), None)
        chebi_id = next((x for x in data_xrefs if x.startswith("CHEBI:")), None)
        hmdb_id = next((x for x in data_xrefs if x.startswith("HMDB")), None)
        return {
            "PubChem": query_id,
            "KEGG": kegg_id,
            "HMDB": hmdb_id,
            "ChEBI": chebi_id,
            "InChI": data_inchi.get("InChI"),
            "InChIKey": data_inchi.get("InChIKey")
        }
    except Exception as e:
        logging.warning(f"PubChem no data for {query_id}: {e}")
        return None

############################ CROSS-REFERENCE INTEGRATION ##############################

def integrate_crossrefs(ids_db_list):
    """Joins results from LipidMaps, PubChem and UniProt for all IDs."""
    # Classify IDs
    lm_ids,uni_ids,pchem_cids = classify_ids(ids_db_list)

    # Extract DataBase/s which will not be mapped
    mapped_db = set(LM) | set(UPROT) | set(PCHEM)
    databases = np.unique(ids_db_list[:, 1])
    missing_db = [db for db in databases if db not in mapped_db and db != 'nan']   
    # Totales por tipo:
    logging.info(
            f"\tLIPID MAPS: {len(lm_ids)} | "
            f"UniProt: {len(uni_ids)} | "
            f"PubChem: {len(pchem_cids)} | "
            f"Databases that will not be mapped: {len(missing_db)} > " 
            f"{', '.join(missing_db) if missing_db else 'N/A'}"
        )
    
    results = {}
    # Query LipidMaps DB (for LM ID)
    for query_id in lm_ids:
        info = fetch_lipidmaps_info(query_id)
        if not info:
            continue
        results[query_id] = info

    # Query UniProt DB
    for query_id in uni_ids:
        if query_id not in results:
            info = fetch_refseq_from_uniprot(query_id)
            if info:
                results[query_id] = info
            else:
                logging.info(f"NO xrefs for {query_id}.")
            # time.sleep(0.2)

    # Query PubChem for uncrossed numeric ids
    for query_id in pchem_cids:
        info = fetch_pubchem_info(query_id)
        if info:
            results[query_id] = info
        else:
            logging.info(f"NO xrefs for {query_id}.")
            # time.sleep(0.3)
            
    # df_results = pd.DataFrame([{"ID": k, **v} for k, v in results.items()])
    df_results = pd.DataFrame.from_dict(results, orient="index")
    df_results.index.name = "ID"
    df_results = df_results.reset_index()
    return df_results


############################ MERGING DataFrame Results ###############################

def merge_df_with_crossrefs(df, df_results):
    """Update the DataFrame with the cross-references found."""
    outer_m_df = pd.merge(df, df_results, on ="ID", how="left")
    return outer_m_df


################################ READ/SAVE FILE #######################################

def read(input_file):
    # Read input + removing empty rows and columns (".dropna")
    df_all = (
        pd.read_csv(input_file, sep='\t', encoding="cp1252")
        .dropna(axis=0, how='all')
        .dropna(axis=1, how='all')
    )
    return df_all

def extract_ids(df):
    # Normalize info in DataBase column -> pass to string and convert to lowercase
    df["DataBase"] = df["DataBase"].astype(str).str.strip().str.lower()

    # Extract ID and DataBase columns
    df_ids_db = df[["ID","DataBase"]]
    ids_db_list = df_ids_db.to_numpy()

    # Check for duplicate IDs
    duplicated = df_ids_db["ID"].duplicated(keep=False)

    repeated_ids = df_ids_db.loc[duplicated, "ID"].dropna().unique() # remove NaN
    if len(repeated_ids) > 0:
        logging.warning(
            f"ALERT: The following IDs are repeated: "
            f"{', '.join(map(str, repeated_ids))}")
        
    return ids_db_list


def save(base_filename, output_folder, final_df):
    # Guardar el df resultante en tsv
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    base_filename_upd = base_filename.replace(".txt", "_updated.txt")
    output_path = os.path.join(output_folder,base_filename_upd)
    final_df.to_csv(output_path, sep='\t')
    return output_path


############################## MAIN CONTROLADOR #######################################

def main(input_file, output_folder):
    total_time = []
    start = time.perf_counter()
    
    base_filename = os.path.basename(input_file)
    logging.info(f" --- Procesando archivo: {base_filename} ---")

    # Read input
    df_all = read(input_file)
    ids_db_list = extract_ids (df_all)

    # Integrate and merge results
    df_results = integrate_crossrefs(ids_db_list)
    final_df = merge_df_with_crossrefs(df_all, df_results)

    # Guardar el df resultante en tsv
    output_path = save(base_filename, output_folder, final_df)

    # Tiempo de ejecución y registro de la ejecución
    end = time.perf_counter()
    total_time.append(end - start)
    logging.info(f"     Tiempo ejecución hoja: {end - start:.2f} s")
    logging.info(f" ** Archivo guardado en: {output_path} **\n")


################################ ENTRY POINT ##########################################

if __name__ == "__main__":
    tsv_files = [f for f in os.listdir(INPUT_FOLDER) 
                 if f.endswith('.txt') or f.endswith('.tsv')]
    for file in tsv_files:
        INPUT_FILE = os.path.join(INPUT_FOLDER, file)
        main(INPUT_FILE, OUTPUT_FOLDER)