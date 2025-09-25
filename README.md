# Pathways-visualization-tool

**WIP (prototype).**

Converts tabular pathway specs (CSV) into **GPML** so they can be opened/saved in **PathVisio/WikiPathways**. Current goal: generate GPML with correct DataNodes (metabolites/enzymes) and Interactions (conversion/catalysis) with basic layout.

> Status: under active development. APIs, file names, and outputs may change.

## What it does (today)
- Parse a CSV describing metabolites/enzymes/interactions.
- Emit a `.gpml` file that opens in PathVisio.
- Include simple positioning & sizes for nodes; node positioning is being refined.

## Repo layout
```
src/ # library / main script
  to_gpml.py # current implementation (work in progress)
examples/
  data/ # example CSV inputs used for tests/demos
  gpml/ # GPML outputs generated from those CSVs
LICENSE
README.md
```
