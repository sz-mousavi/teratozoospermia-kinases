# Data

Raw GEO files are **not** stored in this repository. Download them from GEO and
save them in this folder (`data/`) with exactly these names.

## 1. Expression matrices (one per study)

| Study | Platform | Samples | Save as |
|-------|----------|---------|---------|
| [GSE6967](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6967) | GPL2507 (Illumina Human-6 BeadChip) | 13 (8 teratozoospermic, 5 normozoospermic) | `GSE6967_series_matrix.txt` |
| [GSE6968](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6968) | GPL2700 (Illumina HumanRef-8 BeadChip) | 10 (6 teratozoospermic, 4 normozoospermic) | `GSE6968_series_matrix.txt` |
| [GSE6872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6872) | GPL570 (Affymetrix U133 Plus 2.0) | 21 (8 teratozoospermic, 13 normozoospermic) | `GSE6872_series_matrix.txt` |

How to make each `*_series_matrix.txt`:

1. On the GEO series page, download the *Series Matrix File* (`..._series_matrix.txt.gz`) and unzip it.
2. Keep **only the data table**: delete every line up to and including
   `!series_matrix_table_begin`, and delete the last line `!series_matrix_table_end`.
   The first line of the saved file must be the header
   (`"ID_REF"  "GSM..."  "GSM..."  ...`).

## 2. Platform annotation tables

| Platform | Save as |
|----------|---------|
| [GPL2507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL2507) | `GPL2507.txt` |
| [GPL2700](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL2700) | `GPL2700.txt` |
| [GPL570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570) | `GPL570.txt` |

On each platform page use *Download full table...* and save the tab-delimited
table. Only the `ID` and `GB_ACC` (GenBank/RefSeq accession) columns are used.
The first line of the saved file must be the column header (`ID  ...  GB_ACC ...`),
without the `^`, `#` and `!` header lines.

## Sample order

The scripts label samples by column position (see `gr` in the script / the
"Phenotype labels" section of the notebook). The order of the three matrices is
kept exactly as GEO provides it, so the labels are, in order:

| Columns of merged matrix | Study | Labels |
|---|---|---|
| 1-13 | GSE6967 | 5 terato, 1 normo, 3 terato, 4 normo |
| 14-23 | GSE6968 | 6 terato, 4 normo |
| 24-44 | GSE6872 | 13 normo, 8 terato |

Total: 44 samples (22 teratozoospermic, 22 normozoospermic).

Do not commit the raw files (`data/*` is ignored by `.gitignore`).
