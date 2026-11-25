# Dysregulated kinases in teratozoospermia

Code accompanying:

> Mousavi SZ, Hadizadeh M, Mohammad Soltani B, Totonchi M.  
> Dysregulated kinase expression in teratozoospermia and implications for male infertility: an integrated gene expression study.  
> *Journal of Reproduction & Infertility*. 2025;26(4):224.  
> https://doi.org/10.18502/jri.v26i4.21088

## Method

Meta-analysis of three public microarray datasets comparing teratozoospermic and normozoospermic sperm RNA profiles.

Steps: log2 transformation, probe annotation and merge on GenBank/RefSeq accession (`GB_ACC`), per-gene mean over probes, ComBat batch correction (batch = study), quantile normalization, limma differential expression (teratozoospermia vs normozoospermia; DEGs: |log2FC| > 1 and FDR < 0.01).

## Datasets and sample groups

| Study | Platform | Teratozoospermic | Normozoospermic |
|-------|----------|:---:|:---:|
| [GSE6967](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6967) | GPL2507 | 8 | 5 |
| [GSE6968](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6968) | GPL2700 | 6 | 4 |
| [GSE6872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6872) | GPL570 | 8 | 13 |
| **Total** | | **22** | **22** |

The merged matrix has 44 samples. Group labels are assigned by column order (variable `gr`), matching the order of the samples in the GEO series matrices.

## Data preparation

Raw GEO files are not stored in the repository. Download the three series matrices and the three platform tables from GEO, and save them in `data/` as:

```
data/
  GSE6967_series_matrix.txt   GPL2507.txt
  GSE6968_series_matrix.txt   GPL2700.txt
  GSE6872_series_matrix.txt   GPL570.txt
```

Each series matrix is the GEO *Series Matrix File* reduced to its data table (the `!`-prefixed header lines are removed so the first line is `"ID_REF" "GSM..." ...`). Each platform table comes from *Download full table...* on the GEO platform page; only the `ID` and `GB_ACC` columns are used. Step-by-step details: [`data/README.md`](data/README.md).

## Repository layout

- `meta_analysis_teratozoospermia/teratozoospermia_kinase_meta.Rmd` — annotated notebook of the analysis (canonical, with figures)
- `meta_analysis_teratozoospermia/teratozoospermia_kinase_meta.md` — rendered GitHub version of the notebook
- `scripts/teratozoospermia_kinase_meta.R` — the same analysis as a plain script
- `data/` — dataset download instructions (raw files are not committed)

## Requirements

Developed with R 4.5.0. Package versions used:

| Package | Version |
|---------|---------|
| limma | 3.66.0 |
| sva | 3.58.0 |
| data.table | 1.18.6.1 |
| ggplot2 | 4.0.3 |
| pheatmap | 1.0.13 |
| dplyr | 1.2.1 |
| rmarkdown (notebook only) | — |

Install:

```r
install.packages(c("data.table", "ggplot2", "pheatmap", "dplyr", "rmarkdown"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "sva"))
```

The full `sessionInfo()` is at the end of the rendered notebook.

## How to reproduce

Run everything from the repository root.

1. Prepare `data/` as described above.
2. Run the script:

   ```sh
   Rscript scripts/teratozoospermia_kinase_meta.R
   ```

   or, to re-knit the notebook:

   ```r
   rmarkdown::render(
     "meta_analysis_teratozoospermia/teratozoospermia_kinase_meta.Rmd",
     output_format = "github_document",
     knit_root_dir = getwd()
   )
   ```

Tables and figures are written to `results/` (e.g. `terato-normo.csv`, `Up-DEGs1.csv`, `Down-DEGs1.csv`, `all-DEGs1.csv`, `DEGs-Matrix.csv`). Both routes use the same code, so they give the same DEG lists.

## License

MIT, see [`LICENSE`](LICENSE).

## Citation

Please cite the paper above.
