# saRNA Prediction API

A REST API for predicting candidate small activating RNA (saRNA) sequences from human genes. Based on the hypothesis that messenger RNA fragments can escape NMD-mediated degradation and act as endogenous transcriptional activators on promoters of functionally similar genes.

## Requirements

### R (version 4.3.3 or higher)
```r
install.packages(c("plumber", "httr", "jsonlite", "xml2", 
                   "stringr", "dplyr", "reticulate"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("Biostrings")
```

### Python (version 3.12 or higher)
```bash
pip install pymongo pandas
```

### External tools
- [EMBOSS](http://emboss.open-bio.org/) (version 6.6.0.0)
```bash
sudo apt-get install emboss
```

- [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) (version 2.5.1)
```bash
sudo apt-get install vienna-rna
```

### MongoDB (version 8.2.6 or higher)
```bash
sudo apt-get install mongodb
```

## Database setup

First, download the human GO annotations file from EBI:
```bash
wget https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/25.H_sapiens.goa \
     -P data/
```

Then load it into MongoDB:
```bash
python3 data/load_data.py
```

This processes the GOA file, filters annotations with the `involved_in` 
relationship, removes duplicates and loads the data into the 
`saRNAs.go_annotations` MongoDB collection. The process inserts one 
document per gene with the following structure:

```json
{
  "gene_name": "ATOH1",
  "go_id": ["GO:0061564", "GO:0007417", "GO:0014014", ...]
}
```
## Project structure

```
saRNAs/
├── plumber.R          # Main API
├── run_api.R          # API launcher  
├── functions.R        # Pipeline functions
├── similarity.py      # GO similarity module
├── banner_github.png  # Banner for README
├── README.md
├── data/
│   ├── load_data.py   # MongoDB data loader
│   └── 25.H_sapiens.goa  # GOA file (download separately, **see setup**)
└── promoters/
    └── Hs_EPDnew.bed  # EPD promoter database (GRCh38/Gencode v28)
```
## Usage

### Launch the API
```bash
Rscript run_api.R
```

The API will be available at `http://localhost:8000`.  
Interactive documentation (Swagger): `http://localhost:8000/__docs__/`






