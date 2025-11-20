# Pathways-visualization-tool

Python workflows to support a metabolomics pathway project focused on bile acids and cardiovascular disease.\
The repository currently includes:

- A cross‑referencing pipeline to harmonize metabolite and protein identifiers across major databases (LIPID MAPS, PubChem, UniProt and others).  
- A prototype converter from tabular pathway descriptions (CSV/TSV) into GPML, the XML format used by PathVisio and WikiPathways for pathway visualization and sharing.

These components are being developed as part of a Master’s thesis on the visualization and programmatic querying of metabolomic pathways.

## Project goals

- Make manually reconstructed bile acid pathways machine‑readable and reusable.
- Provide a reproducible workflow to:
  - Clean and annotate metabolite/protein IDs via public REST APIs.
  - Convert curated pathway tables into GPML for visualization in PathVisio and potential upload to WikiPathways.
- Enable later integration with multi‑omics tools (e.g. TurbOmics) for pathway‑level analysis of metabolomics data.

## Repository structure

```
├── src/
│ ├── DB_crosslink.py # Cross‑referencing pipeline for metabolite & protein IDs
│ └── to_gpml.py # CSV → GPML prototype for pathway visualization
├── examples/
│ ├── data/ # Example input tables (metabolites, enzymes, reactions)
│ └── gpml/ # GPML outputs generated with the prototype
├── LICENSE
└── README.md
```

The code is a work in progress; interfaces and file formats may change as the thesis evolves.
