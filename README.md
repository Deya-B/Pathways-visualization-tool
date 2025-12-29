# Pathways Visualization Tool for Metabolomics


Python tools for:
- Harmonizing metabolite/protein IDs across databases.
- Converting pathway ID/relations tables (TSV/CSV) into GPML for visualization and analysis in PathVisio and WikiPathways.

---

## Components

- **Cross‑referencing pipeline** (`DB_crosslink.py`):  
  Enriches metabolite and protein tables with cross‑references from LipidMaps, PubChem, UniProt, KEGG, HMDB, ChEBI, RefSeq, and others via public REST APIs.

- **TSV‑to‑GPML pathway builder** (`src/pathways_tool/` + `src/main.py`):  
  Modular library and CLI that reads ID metadata and relations tables, computes an automatic graph layout, and exports a GPML 2013a pathway file.

The resulting `.gpml` files can be used directly for computational analysis or opened in **PathVisio**, refined manually, and exported for use in **WikiPathways** and other downstream analyses.

For detailed usage, input formats, and API endpoints, see [`src/README.md`](src/README.md).