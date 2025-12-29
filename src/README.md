# Source code overview

This folder contains the main Python pipelines used in the project:

- [`DB_crosslink.py`](#cross‑referencing-pipeline): cross‑referencing pipeline for metabolite and protein identifiers.
- [`to_gpml.py`](#tsv‑to‑gpml-pathway-builder): prototype converter from tabular pathway descriptions to GPML.

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

- Tab‑separated file (`.tsv` / `.txt`), encoded in UTF-8.
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
2. **Run the script**:
     ```
     python DB_crosslink.py -i path/to/input_folder -o path/to/output_folder
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
| `-l` / `--log`     | Logging level (DEBUG/INFO/WARNING/ERROR).<br>Set to INFO by default, DEBUG can be requested with<br>`-l DEBUGG` to obtain more information.|
### Dependencies

- Python 3.9+
- `pandas`, `numpy`, `requests`.

Install via:
```bash
pip install pandas numpy requests
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


## TSV‑to‑GPML pathway builder

The `pathways_tool` package provides a **modular library** and **command‑line interface** to build GPML pathway files from two TSV/CSV tables: one with ID metadata (metabolites, pathways, enzymes) and one with source–target relationships.

- Core modules:
  - `pathways_tool.parser`: builds **Node** and **Interaction** objects from pandas DataFrames.
  - `pathways_tool.layout`: computes a 2D layout (layers, positions, anchors) for nodes and interactions.
  - `pathways_tool.xml_builder`: generates GPML 2013a XML from the in‑memory graph.
  - `pathways_tool.config`: loads YAML configs, validates required fields, resolves paths, and sets up logging.
  - `pathways_tool.cli`: orchestrates reading inputs, building the graph, layout, and writing the GPML file.

### Input files

The GPML builder expects:

- **ID metadata table** (TSV/CSV):
  - Contains identifiers (ID), database names and display labels.
  - Used to create `Node` objects for metabolites, enzymes, and pathway references (i.e. nodes that point to other pathways), with proper labels and xrefs.

- **Relations table** (TSV/CSV):
  - Describes reactions as source → target plus optional catalysts.
  - Used to build:
    - **Conversion** interactions (substrate → product)
    - **Catalysis** interactions (enzyme → reaction anchor)

Example tables and schemas are provided under `examples/data/`.

### YAML configuration

The command‑line entry point uses a YAML config to describe titles, organism, input files, and output file.

Minimal example:

```yaml
# config.yaml
name: "MurineCDCA-MCA_MDCA"

pathway_title: "{name}"
organism: "Homo sapiens"

input:
  folder: "C:/path/to/examples/data/"
  id_data_file: "{name}.tsv"
  relations_file: "{name}_relationships.tsv"
  delimiter: "\t"

output:
  folder: "C:/path/to/examples/data/"
  filename: "{name}.gpml"

logging:
  level: "INFO"
  format: "%(levelname)s: %(message)s"
```

The `name` field can be reused inside strings via `{name}`, and `load_config` expands these templates and checks that required keys are present.

### Command‑line usage

From the repository root:

```bash
python src/main.py --config path/to/config.yaml
```

This will:

- Load and validate the YAML configuration into a `GPMLConfig` dataclass.
- Configure logging according to `logging.level` and `logging.format`.
- Read the ID and relations tables with robust CSV/TSV loading (`io_utils.read_csv`).
- Build nodes and interactions with `Parser` from `pathways_tool.parser`.
- Compute an automatic pathway layout with `Layout` from `pathways_tool.layout`.
- Generate a GPML XML tree with `XMLBuilder` and write it to the configured output file.

The resulting `.gpml` file can be used directly for computational analysis, or it can be opened in **PathVisio**, refined with manual edits to improve visualization, and then exported for upload and use in **WikiPathways** and other downstream computational analyses.

### Library usage from Python

The same functionality can be used programmatically:

```python
from pathways_tool.config import load_config
from pathways_tool.cli import run_from_config

cfg = load_config("path/to/config.yaml")
run_from_config(cfg)
```

Or, if you already have pandas DataFrames:

```python
from pathways_tool.parser import Parser
from pathways_tool.layout import Layout
from pathways_tool.xml_builder import XMLBuilder
import xml.etree.ElementTree as ET

builder = Parser(id_data_df, relations_df)
layout = Layout(builder.nodes, builder.interactions)
layout.layout_positions()
layout.layout_anchors()
layout.layout_catalysis()

xml_builder = XMLBuilder(
    title="My pathway",
    organism="Homo sapiens",
    nodes=builder.nodes,
    interactions=builder.interactions,
)
tree = ET.ElementTree(xml_builder.to_etree())
tree.write("my_pathway.gpml", encoding="utf-8", xml_declaration=True)
```
