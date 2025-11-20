# Source code overview

This folder contains the main Python pipelines used in the project:

- [`DB_crosslink.py`](#cross‑referencing-pipeline): cross‑referencing pipeline for metabolite and protein identifiers.
- [`to_gpml.py`](#gpml-converter-prototype): prototype converter from tabular pathway descriptions to GPML.

---

## Cross‑referencing Pipeline

`DB_crosslink.py` annotates metabolite and protein identifiers from **LipidMaps**, **UniProt**, and **PubChem**, and adds cross‑references to several metabolite databases plus RefSeq accessions for proteins.

### What it does

- Reads one or more tab‑separated files (`.tsv` / `.txt`) from an input folder.
- Extracts `ID` / `DataBase` pairs and classifies them into:
  - LIPID MAPS (metabolites or proteins),
  - UniProt accessions,
  - PubChem CIDs.
- Queries public REST APIs:
  - LIPID MAPS REST service for metabolite and protein records (KEGG, PubChem CID, HMDB, ChEBI, RefSeq, UniProt, InChI, InChIKey).
  - UniProt REST API for RefSeq cross‑references.
  - PubChem PUG REST API for InChI / InChIKey and external IDs (LIPID MAPS, KEGG, HMDB, ChEBI) via Registry IDs.
- Merges the retrieved annotations back into the original table, preserving the input column order.
- Checks for duplicated identifiers across key ID columns (PubChem, HMDB, ChEBI, UniProt) and logs warnings.
- Writes updated files to an output folder with the suffix `_updated.tsv`.

### Input format

- Tab‑separated file (`.tsv` / `.txt`).
<!-- , encoded in cp1252. -->
- Required columns:
  - `ID`:  
    - Metabolites: LIPID MAPS ID or PubChem CID.  
    - Proteins/enzymes: LIPID MAPS protein ID or UniProt accession.  
    - Multiple UniProt IDs can be separated by `;` in one cell (e.g. `P12345;Q8ABC1`).
  - `DataBase`, one of:
    - `LipidMaps`
    - `PubChem`
    - `UniProt`.

- Additional columns (e.g. `Common_Name`, `Synonyms`, comments) are allowed and will be kept untouched.
- To receive specific annotations, include appropriately named empty columns, e.g.:
  - `KEGG`, `PubChem`, `HMDB`, `ChEBI` (xref DBs),
  - `InChIKey`, `InChI`,
  - `RefSeq` (for protein RefSeq accessions).

### Typical workflow

1. **Prepare input tables** with `ID` and `DataBase` plus any extra columns you want to populate.
2. **Configure input/output folders** via:
   - a YAML config file, or
   - CLI arguments if `--config` is not required.
3. **Run the script**:
   - With config file:
     ```
     python DB_crosslink.py -c config.yaml
     ```
   - With manual paths:
     ```
     python DB_crosslink.py \
       -i path/to/input_folder \
       -o path/to/output_folder \
       -v \
       -l INFO
     ```
4. **Inspect logs** for:
   - unmapped databases,
   - duplicated IDs across annotation columns,
   - API timeouts or missing records.

#### CLI Options

| Flag               | Description                              |
| ------------------ | ---------------------------------------- |
| `-i` / `--input`   | Folder containing input TSV/TXT files    |
| `-o` / `--output`  | Output directory for updated files       |
| `-c` / `--config`  | YAML configuration file                  |
| `-v` / `--verbose` | Verbose debugging mode                   |
| `-l` / `--log`     | Logging level (DEBUG/INFO/WARNING/ERROR) |

### Dependencies

- Python 3.9+
- `pandas`, `numpy`, `pyyaml`, `requests`.

Install via:
```bash
pip install pandas numpy pyyaml requests
```

### Input/Output examples:

**Input example**

If our input file `example.txt` is as follows:

| ID | DataBase | Common_Name | Synonyms | KEGG | PubChem | HMDB | ChEBI | InChIKey| Comments |
|---|---|---|---|---|---|---|---|---|---|
| LMST01010001 | LipidMaps | Cholesterol |  |  |  |  |  |  | example 1 |
| LMST04010032 | LipidMaps | Chenodeoxycholic acid | CDCA |  |  |  |  |  | example 2 |
| ... |  |  |  |  |  |  |  |  | example $n$ |

**Output example**

The script generates this `example_updated.txt`:


| ID | DataBase | Common_Name | Synonyms | KEGG | PubChem | HMDB | ChEBI | InChIKey | Comments | 
|---|---|---|---|---|---|---|---|---|---|
| LMST01010001 | LipidMaps | Cholesterol |  | C00187 | 5997 | HMDB0000067 | 16113 | HVYWMOMLDIMFJA-DPAQBDIFSA-N | example 1 |
| LMST04010032 | LipidMaps | Chenodeoxycholic acid | CDCA | C02528 | 10133 | HMDB0000518 | 16755 | RUDATBOHQWOJDD-BSWAIDMHSA-N | example 2 |
| ... |  |  |  |  |  |  |  |  | example $n$  |

### Endpoints used for each API:

| DB        | Endpoint used                                         |
| --------- | ----------------------------------------------------- |
| LipidMaps | `/rest/protein/lmp_id/{id}/all`                       |
|           | `/rest/protein/uniprot_id/{id}/all`                   |
|           | `/rest/compound/lm_id/{id}/all`                       |
| UniProt   | `https://rest.uniprot.org/uniprotkb/{id}?format=json` |
| PubChem   | `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{id}/property/InChI,InChIKey/JSON`|
|           | `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{id}/xrefs/RegistryID/JSON`|


---

## GPML converter prototype

`to_gpml.py` is a **work‑in‑progress** prototype that converts a CSV description of a pathway into a **GPML** file that can be opened and edited in **PathVisio** and potentially uploaded to **WikiPathways**.

### Current capabilities

- Parse a CSV file listing:
  - metabolites and products as nodes,
  - enzymes as separate nodes,
  - reactions / conversions (substrate → product),
  - catalytic relationships (enzyme catalyzes conversion).
- Build an internal pathway graph:
  - `Node` objects for metabolites and enzymes, with labels and database xrefs.
  - `Interaction` objects for:
    - **Conversion** edges (substrate → product),
    - **Catalysis** edges (enzyme → reaction anchor).
- Apply an automatic layout:
  - BFS‑style layering for metabolite/product nodes.
  - Placement of enzymes near the anchor point(s) of the reactions they catalyze.
  - Computation of a GPML “board” size based on node positions.
- Export a GPML 2013a‑style XML file with:
  - `DataNode` elements for metabolites and enzymes.
  - `Interaction` elements for conversions and catalysis, including anchors.
  - Basic node sizes, fonts, and colors.

### Expected input

- A CSV file (default delimiter `\t`; the parser can be configured) containing at least:
  - Node labels (metabolites, products, enzymes).
  - Columns describing reaction source, target, and catalytic enzymes (exact column names depend on the current prototype and may still change).  
- Example input tables are provided under `examples/data/`.

### Usage example

Minimal example (hard‑coded in the prototype at the moment):

```bash
python to_gpml.py
```

This will:

- Read an example CSV such as `examples/data/rutamedia.tsv` (path currently set in the script).
- Build a `Pathway` object.
- Run the internal layout.
- Save a `.gpml` file that can be opened in PathVisio.

Open the resulting GPML in **PathVisio** and refine the layout or annotations there as needed.

### Status and roadmap

- Status: **prototype**, not yet a stable library.
- Known limitations:
  - CSV schema is still evolving.
  - Limited support for compartments, complexes, and advanced GPML features.
- Planned improvements:
  - Stable CLI entry point (input CSV, title, organism, output path).
  - More flexible parsing of column names for metabolites/enzymes/reactions.
  - Better layout heuristics and optional manual overrides.
