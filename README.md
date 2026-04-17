# GDSC2 Data Integration Pipeline

## Overview
This project provides a reproducible pipeline for integrating, cleaning, and standardizing pharmacogenomic data from the GDSC2 dataset.

The objective is to generate a **machine learning–ready dataset** for cancer drug sensitivity prediction.

---

## Key Features
- Integration of drug response, compound annotations, and cell line metadata
- Harmonization of drug targets and biological pathways
- Removal of missing and duplicate observations
- Drug-specific outlier detection
- Missing data analysis
- Correlation analysis of pharmacological metrics
- Analytical validation through visualization

---

## Output
- `outputs/GDSC2_Final_Master_Clean.csv`: curated dataset ready for ML applications

---

## Methods Summary
1. Data integration using COSMIC_ID and DRUG_NAME
2. Annotation harmonization (TARGET / PATHWAY)
3. Filtering missing LN_IC50 values
4. Collapsing replicate measurements
5. Outlier detection (mean ± 3 SD per drug)
6. Correlation structure validation
7. Data quality assessment

---

## Figures
The pipeline generates key validation plots:
- Distribution of LN_IC50
- Correlation between LN_IC50 and AUC
- Missing data profile
- Correlation matrix of response variables

---

## Data Availability

The dataset used in this project is publicly available at:

https://www.cancerrxgene.org/

Download and place the following files in the `data/` folder:

- GDSC2-dataset.csv
- Compounds-annotation.csv
- Cell_Lines_Details.xlsx

---

## Reproducibility

This project was developed and tested using **Google Colab**.

However, it is fully reproducible in any standard R environment.

### Requirements
- R (>= 4.0)
- tidyverse
- readxl
- corrplot

---

## Project Structure

data/  
scripts/  
outputs/  

---

## Future Work
- Feature engineering (multi-omics integration)
- Machine learning modeling
- Biomarker discovery
- Deep learning approaches

---

## Author
David Villafañe
