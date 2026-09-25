
# Replication code for: A map of RBP-miRNA regulatory connections for deciphering the mechanisms of transcriptome remodeling in cancer. 

R project and notebooks required to reproduce the analyses, figures, and tables presented in the publication.

## Repository structure

```text
RBP-miRNA-manuscript/
├── .gitattributes
├── .gitignore
├── .Rprofile
├── LICENSE
├── LICENSE-DATA
├── README.md
├── RBP-miRNA-manuscript.Rproj
├── renv.lock
├── renv/
│   ├── .gitignore
│   ├── activate.R
│   ├── settings.json
│   └── staging/
├── CRISPR_RBP_KO/
│   ├── analysis/
│   └── output/
├── DDX55_miR29_overexpression_inhibition/
│   ├── analysis/
│   └── output/
├── ENCODE_RBP_KD/
│   ├── analysis/
│   └── output/
├── ENCODE_tissues_cell_lines/
│   ├── analysis/
│   └── output/
└── TCGA/
    ├── analysis/
    └── output/
```

Notebooks and their outputs are organized by dataset/experiment. Each `analysis/` 
directory holds the notebooks that reproduce the figures, in both Rmd (source) 
and HTML (rendered) formats. Each `output/` directory holds the generated figures, 
tables, and data objects required for downstream analyses. All notebooks read data from the `data/` directory and write only to their respective `output/` directories.

## System Requirements and Reproducibility

- R Version: 4.4.1
- Package Management: This project uses the `renv` package to ensure computational reproducibility. 

## Installation

1. Download and unzip this repository.
2. Open the `RBP_miRNA_map_manuscript_code.Rproj` file in RStudio.
3. Run `renv::restore()` in the R console to install the exact package versions required for this analysis. Note: If `renv::restore()` fails to build fgsea from source, run `BiocManager::install("fgsea", update = FALSE)` to install it as a binary instead. Confirm with `packageVersion("fgsea”)` that it matches the version in renv.lock, then re-run `renv::status()` to confirm the project is consistent.
4. The `data/` directory does not exist yet. Input data needs to be dowloaded by the user and placed inside `data/` in the project directory to run the notebooks. This can be done in a single step by running the first notebook `ENCODE_RBP_KD/analysis/00_download_input_data.Rmd`. Alternatively, manually download the input data `code_input_data.zip` from Zenodo (https://zenodo.org/records/22679902), unzip the file and place all the files it contains inside a `data/` directory in the project main directory. 

## Running

Open notebooks in RStudio and run the notebooks in the `analysis/` subfolders sequentially. The notebook order is indicated by their names (“01_”, “02_”, etc.). Start with the notebook whose name starts with “01_” (assuming that the input data was already downloaded, see the installation steps above). 

## Figure map

| Notebook | Manuscript Fig. |
|---|---|
| 01_ENCODE_RBP_KD_batch_correction.Rmd | Supplementary Fig. 3 |
| 02_ENCODE_RBP_KD_shRNA_offTargetEffects.Rmd | Supplementary Fig. 5b |
| 03_ENCODE_RBP_KD_miRNA_activity_no_corr_shRNA_offTargetEffects.Rmd | Supplementary Fig. 4, Supplementary Table 2 |
| 04_ENCODE_RBP_KD_miRNA_activity.Rmd | Fig. 2a |
| 05_ENCODE_tissues_cell_lines_stability_miRNA_activity.Rmd | |
| 06_ENCODE_tissues_cell_lines_corr_RBP_expr_miR_act.Rmd | Fig. 2c |
| 07_CRISPR_RBP_KO.Rmd | Fig. 2d, Supplementary Fig. 6d-e, Supplementary Table 3 |
| 08_TCGA_mRNA_stability_bias_correction.Rmd | |
| 09_TCGA_miRNA_activity.Rmd | Fig. 3a, Supplementary Fig. 7, Supplementary Table 4 |
| 10_TCGA_miRNA_activity_vs_abundance.Rmd | Fig. 3b, Supplementary Fig. 8 |
| 11_TCGA_miRNA_transcription_vs_activity.Rmd | Fig. 3c, Supplementary Fig. 9a |
| 12_TCGA_miRNA_transcription_vs_abundance.Rmd | Fig. 3d, Supplementary Fig. 9b |
| 13_TCGA_RBP_expression_vs_miRNA_activity.Rmd | Fig. 4a-b, Fig. 5b, Supplementary Fig. 10b and 12b |
| 14_TCGA_miRNA_activity_vs_purity.Rmd | Supplementary Fig. 10a |
| 15_TCGA_RBPresponsive_miRtargets_enrich_cancer.Rmd | Fig. 4c-e, Supplementary Fig. 11, Supplementary Table 5 |
| 16_TCGA_RBP_expression.Rmd | Fig. 5a, Supplementary Fig. 12a |
| 17_TCGA_DDX55_expression_survival.Rmd | Fig. 5c-d, Supplementary Fig. 10c |
| 18_DDX55_KD_miR-29_mimic_inhibitor.Rmd | Fig. 5e |

Notebooks 05 and 08 have no corresponding figure; they produce intermediate data that feeds the downstream analyses.

## Licence

- Code (notebooks): MIT, see LICENSE.
- Data and figures: CC-BY-4.0, see LICENSE-DATA.

Please cite the paper. 
