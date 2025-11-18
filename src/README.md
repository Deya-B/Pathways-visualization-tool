# Cross‑referencing Pipeline (LipidMaps / PubChem / UniProt)

This repository contains a Python script (`DB_crosslink.py`) that annotates metabolite and protein identifiers (IDs) from **LipidMaps**, **UniProt**, and **PubChem**.It reads TSV files, queries public REST APIs, merges the retrieved annotations, and outputs updated files.

The annotations retrieved for **metabolites** are the *InChI* and *InChIKey* together with the cross-references of IDs with other important metabolite databases: *KEGG*, *PubChem*, *HMDB* and *ChEBI*.\
On the other hand for proteins the RefSeq accession number and UniProt ID are obtained.\

Pipeline steps:
1. Read input files.
2. Extract ID/DataBase pairs.
3. Classify IDs → LipidMaps / UniProt / PubChem.
4. Query APIs.
5. Merge annotations.
6. Check duplicated cross-references.
7. Write output file.

## 1. Requirements

Dependencies required: pandas, numpy, pyyaml and requests

```bash
pip install pandas numpy pyyaml requests
```

## 2. File Structure

Input files must be tab‑separated (`.txt` or `.tsv`) and contain at least:
- ID: Annotation can be performed with IDs from LipidMaps or PubChem in the case of metabolites and LipidMaps or UniProt in the case of proteins/enzymes. Multiple UniProt IDs can be provided in a single cell separated by ";". `Example: P12345;Q8ABC1`. The pipeline will retrieve RefSeq accessions for each ID independently and concatenate the results using ";" between groups and "," inside each group.
- DataBase for the given ID. Accepted DataBases in column:
  - LipidMaps
  - PubChem
  - UniProt

| ID           | DataBase  | xref-DB1 |  xref-DB2 |  ... | InChIKey| InChI|
| ------------ | --------- | --- | --- | --- | --- | --- |
| LMFA01010001 | LipidMaps |||...|||
| P12345       | UniProt   |||...|||
| 1234         | PubChem   |||...|||

- No more data needs to be entered as input. However, it is possible to enter as many columns of information as you require for your data (such as Systematic Name, Synonyms...).
- It is **important** however, that you provide the **name to the columns** for the annotations you wish to obtain.
- The annotations that can be included are:
  - Databases from which to obtain cross-references (xref-DBs) include:
    - LipidMaps, KEGG, PubChem CID, HMDB, ChEBI for metabolites 
    - RefSeq accessions for enzymes
  - InChIKey and InChI

*Example below.*

## 3. Running the Script

### 3.1 Run using a configuration file: 

#### YAML Configuration Example
By creating a simple `config.yaml` with the following information:
```yaml
input_folder: ./relative/input/folder/path
output_folder: c:/User/absolute/output/folder/path
verbose: false
loglevel: "INFO"
```
Then open the terminal, go to the folder where the script and `.yaml` is found, and run:
```bash
python DB_crosslink.py -c config.yaml
```

### 3.2 Run by providing paths manually:

```bash
python DB_crosslink.py \
  -i path/to/input_folder \
  -o path/to/output_folder \
  -v \
  -l INFO
```

#### CLI Options

| Flag               | Description                              |
| ------------------ | ---------------------------------------- |
| `-i` / `--input`   | Folder containing input TSV/TXT files    |
| `-o` / `--output`  | Output directory for updated files       |
| `-c` / `--config`  | YAML configuration file                  |
| `-v` / `--verbose` | Verbose debugging mode                   |
| `-l` / `--log`     | Logging level (DEBUG/INFO/WARNING/ERROR) |


## 4. Input / Output

### 4.1 Input example

If our input file `example.txt` is as follows:

| ID | DataBase | Common_Name | Synonyms | KEGG | PubChem | HMDB | ChEBI | InChIKey | InChI | Comments |
|---|---|---|---|---|---|---|---|---|---|---|
| LMST01010001 | LipidMaps | Cholesterol |  |  |  |  |  |  |  |  | 
| LMST04010032 | LipidMaps | Chenodeoxycholic acid | CDCA |  |  |  |  |  |  |  | 
| ... |  |  |  |  |  |  |  |  |  |  |

### 4.2 Output example

The script generates this `example_updated.txt`:

| ID | DataBase | Common_Name | Synonyms | KEGG | PubChem | HMDB | ChEBI | InChIKey | Comments | InChIKey | InChI | Comments |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| LMST01010001 | LipidMaps | Cholesterol |  | C00187 | 5997 | HMDB0000067 | 16113 | HVYWMOMLDIMFJA-DPAQBDIFSA-N |  |  |  |  |
| LMST04010032 | LipidMaps | Chenodeoxycholic acid | CDCA | C02528 | 10133 | HMDB0000518 | 16755 | RUDATBOHQWOJDD-BSWAIDMHSA-N |  |  |  |  |
| ... |  |  |  |  |  |  |  |  |  |  |  |  |


## 5. Notes

* The pipeline checks for duplicated identifiers across key annotation columns (PubChem, HMDB, ChEBI, UniProt). When duplicates occur, a warning is printed.
* UniProt RefSeq version numbers (`.1`, `.2`, …) are removed.
* Multi‑entries for UniProt IDs written as `ID1;ID2` are processed individually and rejoined.
* Note: Large batches may trigger API rate-limits. The script does not implement automatic backoff. If needed, introduce delays in the API section.

## 6. Endpoints used for each API:

| DB        | Endpoint used                                         |
| --------- | ----------------------------------------------------- |
| LipidMaps | `/rest/protein/lmp_id/{id}/all`                       |
|           | `/rest/protein/uniprot_id/{id}/all`                   |
|           | `/rest/compound/lm_id/{id}/all`                       |
| UniProt   | `https://rest.uniprot.org/uniprotkb/{id}?format=json` |
| PubChem   | `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{id}/property/InChI,InChIKey/JSON`|
|           | `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{id}/xrefs/RegistryID/JSON`|


## 7. License

MIT License.
