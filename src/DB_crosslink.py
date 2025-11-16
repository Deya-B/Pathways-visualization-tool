##################################################################################
#          Cross-referencing pipeline (LipidMaps / PubChem / UniProt)            # 
##################################################################################
# This script:
    # 1. Reads TSV files with "ID" and "DataBase" columns.
    # 2. Sorts the IDs according to their source: LipidMaps, UniProt, or PubChem.
    # 3. Queryes the REST APIs of:
        # • LipidMaps:
                # KEGG, PubChem CID, HMDB, ChEBI, RefSeq_Id, UniProt, InChI, InChIKey
        # • UniProt:
                # Cross-refs RefSeq 
                    # multiple accesses per UID; 
                    # .1/.2 suffixes are removed
        # • PubChem: 
                # InChI, InChIKey, and cross-refs KEGG, HMDB, and ChEBI
    # 4. Integrates the results into a single dataframe.
    # 5. Merges the original dataframe with the new annotation columns.
    # 6. Export the updated file with the suffix "_updated.txt".
##################################################################################
# Pipeline Architecture:
# ├── fetch_lipidmaps_info()     # LipidMaps Query: metabolites and proteins
# ├── fetch_refseq_from_uniprot()# UniProt Query: clean retrieval of RefSeq accessions
# ├── fetch_pubchem_info()       # PubChem Query: InChI/InChIKey + cross-refs KEGG/HMDB/ChEBI
# ├── classify_ids()             # Filter IDs by the database declared in the column
# ├── integrate_crossrefs()      # Merge results from the three sources into a single dict/DF
# ├── merge_df_with_crossrefs()  # Combine the original TSV with the new annotations
# ├── read(), save()             # File input/output
# └── main()                     # General file processing controller
##################################################################################

import os
import pandas as pd
import numpy as np
import requests
import time
import logging
import re

############################# CONFIGURATION AND CONSTANTS ########################
# Folders
INPUT_FOLDER = "C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/pathways_raw"
OUTPUT_FOLDER = "C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/pathways_updated"
# INPUT_FOLDER = "c:/Users/dborrotoa/Desktop/TFM/pathways_raw/"
# OUTPUT_FOLDER = "c:/Users/dborrotoa/Desktop/TFM/pathways_updated/"

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

# UniProt ID Pattern (verified according to official rules)
     # Ref.: https://www.uniprot.org/help/accession_numbers
UNIPROT_PATTERN = (r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z]"
                    "[0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")
   

############################# FETCHERS (APIs) ####################################

def fetch_lipidmaps_info(query_id):
    """Search for information using LipidMaps/UniProt IDs in the LipidMaps 
    REST API."""
    # Obtain context and identifier type segment for DataBase endpoint
        # https://www.lipidmaps.org/rest/
                    # {context}/{input_item}/ >> context and identifier type
                    # {input_value}/{output_item}
    if query_id.startswith("LMP"):      
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
        
        # PROTEIN case: multiple rows in response(r) > use only Row1
        if isinstance(data, dict) and any(k.startswith("Row") for k in data):
            first_row = data.get("Row1") or {}
            return {
                "RefSeq_Id": first_row.get("refseq_id"),
                "UniProt": first_row.get("uniprot_id")
            }
        # METABOLITE/PROTEIN case: flat dict in response
        elif isinstance(data, dict):
            return {
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
    uniprot_ids = [x.strip() for x in str(query_id).split(";")] # rows with ID1;ID2
    for uid in uniprot_ids:
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{uid}?format=json"
            r = requests.get(url, timeout=10)
            r.raise_for_status()
            data = r.json()

            # Extract RefSeq IDs
            refseq_ids = []
            for xref in data.get("uniProtKBCrossReferences", []):
                if xref.get("database") == "RefSeq":
                    refseq_ids.append(xref.get("id"))
            return {
                "UniProt": uid,
                "RefSeq_Id": ",".join(refseq_ids) if refseq_ids else "NaN"
            }
        except Exception as e:
            logging.warning(f"UniProt no data for {query_id}: {e}")
            return None
        

def fetch_pubchem_info(query_id):
    """Search for information by CID in the PubChem API.""" 
    try:    
        # Extract InChI and InChIKey
        url_inchi = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
                     f"{query_id}/property/InChI,InChIKey/JSON")
        r_inchi = requests.get(url_inchi, timeout=10)
        r_inchi.raise_for_status()
        data_inchi = (r_inchi.json()
                      .get("PropertyTable", {})
                      .get("Properties", [{}])[0])
        # Extract crossreferences (xrefs)
        url_xrefs = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
                     f"{query_id}/xrefs/RegistryID/JSON")
        r_xrefs = requests.get(url_xrefs, timeout=10)
        r_xrefs.raise_for_status()
        data_xrefs = (r_xrefs.json()
                        .get("InformationList", {})
                        .get("Information", [{}])[0]
                        .get("RegistryID", []))
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


############################ CROSS-REFERENCE INTEGRATION #########################

def integrate_crossrefs(lm_ids,uni_ids,pchem_cids):
    """Joins results from LipidMaps, PubChem and UniProt for all IDs."""
    results = {}
    # Query LipidMaps DB (for LM ID)
    for query_id in lm_ids:
        info = fetch_lipidmaps_info(query_id)
        if not info:
            continue
        results[query_id] = info

    # Query UniProt DB
    for query_id in uni_ids:
        uids = (query_id).split(";")
        refseq_results = []

        for uid in uids:
            info = fetch_refseq_from_uniprot(uid)
            refseq_raw = info.get("RefSeq_Id") if info else "NaN"
            if refseq_raw == "NaN":
                refseq_results.append("NaN")
                continue
            # Obtain RefSeq parts and remove sufixes
            refseq_raw_parts = [p.strip() for p in refseq_raw.split(",") if p.strip()]
            clean = [p.split(".")[0] for p in refseq_raw_parts]
            refseq_results.append(",".join(clean))
        # Save results
        results[query_id] = {
            "UniProt": query_id,
            "RefSeq_Id": ";".join(refseq_results)
        }
    
    # Query PubChem for uncrossed numeric ids
    for query_id in pchem_cids:
        info = fetch_pubchem_info(query_id)
        if info:
            results[query_id] = info
        else:
            logging.info(f"NO xrefs for {query_id}.")
            # time.sleep(0.3)
    
    df_results = pd.DataFrame.from_dict(results, orient="index")
    df_results.index.name = "ID"
    df_results = df_results.reset_index()
    return df_results


############################# DF INFO EXTRACTION/MERGING #########################

def extract_header (df_all):
    header_list = []
    for col in df_all.columns:
        if not col.startswith("Unnamed"): # remove Unnamed
            header_list.append(col)
    return header_list


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
            f"\n [ALERT]: The following IDs are duplicated:\n"
            f"{','.join(map(str, repeated_ids))}\n")
    return ids_db_list


def classify_ids(ids_db_list):
    """Filtering the "ID" column per DataBase."""
    # Extract ID"s (column 0) where DataBase is matching
    lm_ids = ids_db_list[np.isin(ids_db_list[:, 1], LM), 0]
    uni_ids =  ids_db_list[np.isin(ids_db_list[:, 1], UPROT), 0]
    pchem_cids = ids_db_list[np.isin(ids_db_list[:, 1], PCHEM), 0]
    return lm_ids,uni_ids,pchem_cids


def merge_df_with_crossrefs(df, df_results, col_order):
    """Merge left, join original and API columns, and reorder columns."""
    m_df = pd.merge(df, df_results, on ="ID", how="left")
    m_df = m_df.rename(
        columns={f"{col}_y": col for col in col_order 
                        if f"{col}_y" in m_df.columns}
                        )
    # Order according to headers in the original tsv file
    outer_m_df = m_df[col_order]
    return outer_m_df


################################ READ/SAVE FILE ##################################

def read(input_file):
    # Read input + removing empty rows (".dropna")
    df_all = (
        pd.read_csv(input_file, sep="\t", encoding="cp1252")
        .dropna(axis=0, how="all")
        )
    return df_all


def save(base_filename, output_folder, final_df):
    # Guardar el df resultante en tsv
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    base_filename_upd = base_filename.replace(".txt", "_updated.txt")
    output_path = os.path.join(output_folder,base_filename_upd)
    final_df.to_csv(output_path, sep="\t", index=False)
    return output_path


############################## MAIN CONTROLADOR ##################################

def main(input_file, output_folder):
    total_time = []
    start = time.perf_counter()
    
    base_filename = os.path.basename(input_file)
    logging.info(f" --- Procesando archivo: {base_filename} ---")

    # Read input
    df_all = read(input_file)
    col_order = extract_header(df_all)
    ids_db_list = extract_ids (df_all)

    # Classify IDs
    lm_ids,uni_ids,pchem_cids = classify_ids(ids_db_list)

    # Extract DataBase totals and those which will not be mapped
    mapped_db = set(LM) | set(UPROT) | set(PCHEM)
    databases = np.unique(ids_db_list[:, 1])
    logging.info(
        f" LIPID MAPS: {len(lm_ids)} |"
        f" UniProt: {len(uni_ids)} |"
        f" PubChem: {len(pchem_cids)}"
    )
    missing_db = [db for db in databases if db not in mapped_db and db != "nan"]   
    if missing_db:
        logging.warning(
            f" [ALERT] Databases that will not be mapped: {len(missing_db)}" 
            f" > {','.join(missing_db)}"
        )

    # Integrate and merge results
    df_results = integrate_crossrefs(lm_ids,uni_ids,pchem_cids)
    final_df = merge_df_with_crossrefs(df_all, df_results, col_order)

    # Guardar el df resultante en tsv
    output_path = save(base_filename, output_folder, final_df)

    # Tiempo de ejecución y registro de la ejecución
    end = time.perf_counter()
    total_time.append(end - start)
    logging.info(f" Tiempo ejecución hoja: {end - start:.2f} s")
    logging.info(f" Archivo guardado en:\n\t{output_folder}\n")


################################ ENTRY POINT #####################################

if __name__ == "__main__":
    tsv_files = [f for f in os.listdir(INPUT_FOLDER) 
                 if f.endswith(".txt") or f.endswith(".tsv")]
    for file in tsv_files:
        INPUT_FILE = os.path.join(INPUT_FOLDER, file)
        main(INPUT_FILE, OUTPUT_FOLDER)