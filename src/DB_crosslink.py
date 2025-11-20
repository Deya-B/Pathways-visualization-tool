# python DB_crosslink.py -c config.yaml
###############################################################################
#        Cross-referencing pipeline (LipidMaps / PubChem / UniProt)           # 
###############################################################################
# This script:
    # 1. Reads TSV files with "ID" and "DataBase" columns.
    # 2. Sorts the IDs according to their source:   
    #    LipidMaps, UniProt, or PubChem.
    # 3. Queryes the REST APIs of:
        # > LipidMaps: KEGG, PubChem CID, HMDB, ChEBI, RefSeq_Id, 
        #              UniProt, InChI, InChIKey
        # > UniProt: Cross-refs with RefSeq accession number
        #   - Taking multiple accesses per UID as [ID1;ID2],; which are 
        #     returned as [ID1_RefSeq_ID1, ID1_RefSeq_ID2; 
        #                  ID2_RefSeq_ID1...]
        #   - The version numbers appended to the accession numbers with 
        #     a period (.1/.2...) are removed
        # > PubChem: InChI, InChIKey, and cross-refs KEGG, HMDB, and ChEBI
    # 4. Integrates the results into a single dataframe.
    # 5. Merges the original dataframe with the new annotation columns.
    # 6. Export the updated file with the suffix "_updated.txt".
###############################################################################
# Pipeline Architecture:
# ├── fetch_lipidmaps_info()    # LipidMaps Query: metabolites and proteins         
# ├── fetch_refseq_from_uniprot() # UniProt Query: retrieval of RefSeq
# ├── fetch_pubchem_info()      # PubChem Query: metabolites
# ├── classify_ids()            # Filter IDs by database
# ├── integrate_crossrefs()     # Combine results from all three sources
# ├── merge_df_with_crossrefs() # Add new annotations to original TSV
# ├── read(), save() 
# └── main()                    # General file processing controller
###############################################################################

import os
import pandas as pd
import numpy as np
import requests
import time
import logging
import re
import argparse

########################## CONFIGURATION AND CONSTANTS ########################
parser = argparse.ArgumentParser(description="Process pipeline arguments.")
parser.add_argument("-i", "--input", help="Input folder with TSV files")
parser.add_argument("-o", "--output", help="Output folder path")
parser.add_argument('-l', '--log', 
                    choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'], 
                    default='INFO', help="Log level")
# Parse command line 
args = parser.parse_args()

input_folder = args.input
output_folder = args.output
loglevel = args.log

# Logging Setup
logging.basicConfig(
    level=getattr(logging, loglevel), 
    format="%(levelname)s:%(message)s")

# Global variables with accepted variants
PCHEM = ["pubchem cid", "pubchem"]
LM = ["lipidmaps"]
UPROT = ["uniprot"]
SUPPORTED_DBS = PCHEM+UPROT+['hmbd', 'chebi']
# UniProt ID Pattern (verified according to official rules)
     # Ref.: https://www.uniprot.org/help/accession_numbers
UNIPROT_PATTERN = (r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z]"
                    "[0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")
   
############################ FETCHERS (APIs) ##################################
def safe_get(d, key):
    v = d.get(key)
    return v if isinstance(v, str) else None

def fetch_lipidmaps_info(query_id):
    """Retrieve annotation records from the LipidMaps REST API.

    Parameters
    ----------
    query_id : str
        LipidMaps metabolite ID (LM*), LipidMaps protein ID (LMP*),
        or UniProt accession used as API query key.

    Returns
    -------
    info : dict
        Dictionary containing KEGG, PubChem CID, HMDB, ChEBI, RefSeq,
        UniProt, InChI, and InChIKey fields. Returns None when no data
        are available or the API response is empty.
    """
    # Obtain "context and identifier" type segment for DataBase endpoint
        # https://www.lipidmaps.org/rest/
                    # {context}/{input_item}/ >> "context and identifier"
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
    """Retrieve RefSeq cross-references from the UniProt REST API.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession used to query cross-referenced RefSeq identifiers.

    Returns
    -------
    refseq_ids : list
        List of RefSeq accessions mapped to the UniProt entry.
        Returns an empty list when no RefSeq cross-references exist.
    """
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
    """Retrieve annotation fields from PubChem using a compound CID.

    Parameters
    ----------
    cid : str
        PubChem Compound Identifier used as the query key.

    Returns
    -------
    info : dict
        Dictionary containing LipidMaps, KEGG, HMDB, ChEBI, InChI, and 
        InChIKey extracted from PubChem. Returns None when the CID is 
        invalid or the API response provides no useful data.
    """
    try:    
        # server = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        # Extract InChI and InChIKey
        url_inchi = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
                f"{query_id}/property/InChI,InChIKey/JSON")
        r_inchi = requests.get(url_inchi, timeout=10)
        r_inchi.raise_for_status()
        data_inchi = (r_inchi.json()
                      .get("PropertyTable", {})
                      .get("Properties", [{}])[0])
        # Extract crossreferences (xrefs)
        url_xrefs = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
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
        lm_id = next((x for x in data_xrefs if x.startswith("LM")), None)
        return {
            "PubChem": query_id,
            "LipidMaps": lm_id,
            "KEGG": kegg_id,
            "HMDB": hmdb_id,
            "ChEBI": chebi_id.split(":")[1] if chebi_id else None,
            "InChI": data_inchi.get("InChI"),
            "InChIKey": data_inchi.get("InChIKey")
        }
    except Exception as e:
        logging.warning(f"PubChem no data for {query_id}: {e}")
        return None


########################## CROSS-REFERENCE INTEGRATION ########################

def integrate_crossrefs(lm_ids, uni_ids, pchem_cids):
    """Merge cross-referenced annotations from LipidMaps, UniProt, and PubChem.

    Parameters
    ----------
    lm_ids : list
        IDs matching LipidMaps metabolite or protein patterns.
    uni_ids : list
        IDs matching UniProt accession patterns.
    pchem_cids : list
        IDs recognized as PubChem CIDs.

    Returns
    -------
    df_results : pandas.DataFrame
        DataFrame mapping each input ID to its combined annotation fields.
        When multiple protein IDs are taken as UniProt input: [ID1;ID2],
        these are returned as [ID1_RefSeq_ID1, ID1_RefSeq_ID2;
                               ID2_RefSeq_ID1...].
        Missing values are represented as NaN.
    """
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


########################## DF INFO EXTRACTION/MERGING #########################

def extract_header (df_all):
    """Extract non-empty column names from the input DataFrame.

    Parameters
    ----------
    df_all : pandas.DataFrame
        Full DataFrame loaded from the input TSV file.

    Returns
    -------
    header_list : list
        List of column names excluding unnamed placeholder columns.
    """
    header_list = []
    for col in df_all.columns:
        if not col.startswith("Unnamed"): # remove Unnamed
            header_list.append(col)
    return header_list


def extract_ids(df):
    """Extract ID-DataBase pairs.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing at least 'ID' and 'DataBase' columns.

    Returns
    -------
    ids_db_list : numpy.ndarray
        Array of shape (n, 2) containing ID and DataBase values.
    """
    # Normalize info in DataBase column -> pass to string and convert to lowercase
    df["DataBase"] = df["DataBase"].astype(str).str.strip().str.lower()

    # Extract ID and DataBase columns
    df_ids_db = df[["ID","DataBase"]]
    ids_db_list = df_ids_db.to_numpy()
    return ids_db_list


def classify_ids(ids_db_list):
    """Classify identifiers into LipidMaps, UniProt, or PubChem categories.

    Parameters
    ----------
    id_db_list : list
        List of input identifiers to be classified.

    Returns
    -------
    lm_ids : list
        IDs matching LipidMaps metabolite or protein patterns.
    uni_ids : list
        IDs matching UniProt accession patterns.
    pchem_cids : list
        IDs recognized as PubChem CIDs.
    """
    # Extract ID"s (column 0) where DataBase is matching
    lm_ids = ids_db_list[np.isin(ids_db_list[:, 1], LM), 0]
    uni_ids =  ids_db_list[np.isin(ids_db_list[:, 1], UPROT), 0]
    pchem_cids = ids_db_list[np.isin(ids_db_list[:, 1], PCHEM), 0]
    return lm_ids,uni_ids,pchem_cids


def merge_df_with_crossrefs(df, df_results, col_order):
    """Join the original DataFrame with external annotation fields.

    Parameters
    ----------
    df : pandas.DataFrame
        Input table loaded from TSV file.
    df_results : pandas.DataFrame
        Mapping from ID to extracted annotation records.
    col_order: list
        Header of the input TSV file.

    Returns
    -------
    merged_df : pandas.DataFrame
        DataFrame containing the original columns order with annotated 
        fields.
    """
    m_df = pd.merge(df, df_results, on ="ID", how="left")
    m_df = m_df.rename(
        columns={f"{col}_y": col for col in col_order 
                        if f"{col}_y" in m_df.columns}
                        )
    # Order according to headers in the original tsv file
    output_m_df = m_df[col_order]
    return output_m_df


def duplicated_ids_check(df):
    """Detect duplicates across database ID columns (PubChem, HMDB, ChEBI).
    Duplicated IDs are logged as warnings.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame with cross-referenced data.
    """
    # Identify relevant DB columns from header
    databases = [col for col in df.columns if col.lower() in SUPPORTED_DBS]
    # keep rows where at least one DB column has a value
    df_valid = df[df[databases].notna().any(axis=1)]
    # Check duplicates across those columns
    mask = df_valid.duplicated(subset=databases, keep=False)
    dup_df = df_valid[mask]
    # Extract duplicated IDs grouped per duplicate blocks
    grouped = (
        dup_df.groupby(databases)["ID"]
        .apply(list)              # list of IDs in each duplicate cluster
        .apply(lambda ids: "-".join(map(str, ids)))  # collapse IDs
    )
    # log output
    for merged_ids in grouped:
        logging.warning(
            f"\n\t[ALERT]: Duplicated ID's present: [{merged_ids}]\n"
        )


############################# READ/SAVE FILE ##################################

def read(input_file):
    """Load the input TSV file and remove empty rows.

    Parameters
    ----------
    input_file : str
        Path to the TSV file to be read.

    Returns
    -------
    df_all : pandas.DataFrame
        Parsed DataFrame with empty rows removed.
    """
    df_all = (
        pd.read_csv(input_file, sep="\t", encoding="utf_8")
        .dropna(axis=0, how="all")
        )
    return df_all


def save(base_filename, output_folder, final_df):
    """Write the final DataFrame to a TSV file in the output folder.

    Parameters
    ----------
    base_filename : str
        Original input filename, used to build the output name.
    output_folder : str
        Directory where the updated TSV file will be saved.
    final_df : pandas.DataFrame
        Merged dataFrame containing cross-referenced annotations to write.

    Returns
    -------
    output_path : str
        Full path to the saved TSV file.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    base_filename_upd = base_filename.replace(".tsv", "_updated.tsv")
    output_path = os.path.join(output_folder,base_filename_upd)
    final_df.to_csv(output_path, sep="\t", index=False)
    return None


########################### MAIN CONTROLADOR ##################################

def main(input_file, output_folder):
    """Execute the cross-referencing workflow. Which consists in loading 
    the input TSV, extracting cross-references, merging annotations,
    and writing the output file.

    Parameters
    ----------
    input_file : str
        Path to the TSV file to be read.
    output_folder: str
        Directory where the updated TSV file will be saved.
    """
    total_time = []
    start = time.perf_counter()
    
    base_filename = os.path.basename(input_file)
    logging.info(f" --- Procesando archivo: {base_filename} ---")

    # Read input
    df_all = read(input_file)
    header = extract_header(df_all)
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
            f"\n\t[ALERT]: Databases that will not be mapped: {len(missing_db)}" 
            f" > {','.join(missing_db)}\n"
        )

    # Integrate and merge results
    df_results = integrate_crossrefs(lm_ids, uni_ids, pchem_cids)
    final_df = merge_df_with_crossrefs(df_all, df_results, header)

    # Check duplicated ids + Saving final df to tsv
    duplicated_ids_check(final_df)
    save(base_filename, output_folder, final_df)

    # Tiempo de ejecución y registro de la ejecución
    end = time.perf_counter()
    total_time.append(end - start)
    logging.info(f" Tiempo ejecución hoja: {end - start:.2f} s")
    logging.info(f" Archivo guardado en:\n\t{output_folder}\n")


############################# ENTRY POINT #####################################

if __name__ == "__main__":
    INPUT_FOLDER = input_folder
    OUTPUT_FOLDER = output_folder
    tsv_files = [f for f in os.listdir(INPUT_FOLDER) 
                 if f.endswith(".tsv")]
    for file in tsv_files:
        INPUT_FILE = os.path.join(INPUT_FOLDER, file)
        main(INPUT_FILE, OUTPUT_FOLDER)