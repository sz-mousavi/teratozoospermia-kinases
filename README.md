# Dysregulated kinases in teratozoospermia

Code accompanying:

> Mousavi SZ, Hadizadeh M, Mohammad Soltani B, Totonchi M.  
> Dysregulated kinase expression in teratozoospermia and implications for male infertility: an integrated gene expression study.  
> *Journal of Reproduction & Infertility*. 2025;26(4):224.  
> https://doi.org/10.18502/jri.v26i4.21088

## Method

Meta-analysis of public microarray datasets comparing teratozoospermia and normozoospermia.  
Steps: annotation merge, ComBat batch correction, quantile normalization, limma differential expression.

## Datasets

- [GSE6967](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6967)
- [GSE6968](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6968)
- [GSE6872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6872)

See `data/README.md`.

## Repository layout

- `scripts/teratozoospermia_kinase_meta.R` — analysis script
- `data/` — dataset download notes

## How to run

```r
# from the repository root, after GEO files are in data/
source("scripts/teratozoospermia_kinase_meta.R")

## Data policy

Public datasets only.

## Citation

Please cite the paper above.
