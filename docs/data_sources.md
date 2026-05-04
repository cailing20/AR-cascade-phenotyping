# Data Sources

This repository does not track large primary data files. Use the notes below to
recreate the input layout.

## GWAS Catalog

Current manuscript working file:

```text
/project/CRI/DeBerardinis_lab/lcai/CommonData/GWAScatelog/2026/gwas_catalog_v1.0.2-associations_e115_r2026-04-27_full.tsv
```

Expected repo location:

```text
data/raw/gwas_catalog_v1.0.2-associations_e115_r2026-04-27_full.tsv
```

Used columns:

```text
CHR_ID
CHR_POS
RISK ALLELE FREQUENCY
PVALUE_MLOG
DISEASE/TRAIT
```

Notes:

- Significant associations are filtered at `PVALUE_MLOG > 8`.
- Author-reported risk allele frequency is converted to minor allele frequency
  with `pmin(risk_af, 1 - risk_af)`.
- Records without simple single-chromosome coordinates are excluded.

## ClinVar

Current manuscript working file:

```text
/project/CRI/DeBerardinis_lab/lcai/CommonData/clinvar/May2026/variant_summary.txt.gz
```

Expected repo location:

```text
data/raw/variant_summary.txt.gz
```

Used rows and columns:

- `Assembly == "GRCh38"` for gnomAD matching and allele-frequency plots.
- `Assembly == "GRCh37"` was used only in the older PheWeb overlap workflow,
  which is not part of the cleaned Fig1 rebuild.
- Required columns include `Chromosome`, `Start`, `ReferenceAllele`,
  `AlternateAllele`, `ClinicalSignificance`, `GeneSymbol`, `Name`, and
  `PhenotypeList`.

## gnomAD v4

The extractor streams gnomAD v4 genome VCFs from UCSC:

```text
https://hgdownload.soe.ucsc.edu/gbdb/hg38/gnomAD/v4/genomes
```

Expected repo location for extracted AF files:

```text
data/interim/gnomad_v4_AF/gnomad_v4_AF_chr1.txt.gz
...
data/interim/gnomad_v4_AF/gnomad_v4_AF_chrY.txt.gz
```

Each extracted file has no header and contains:

```text
CHROM  POS  REF  ALT  AF  AF_GRPMAX
```

The current cleaned figure script defaults to plotting `AF`, matching the
original v3 R script. Use `--clinvar-af-column AF_GRPMAX` if the intended panel
should use group-maximum allele frequency.

## PheWeb

The older script referenced lab-local PheWeb TSV files:

```text
/project/CRI/DeBerardinis_lab/lcai/CommonData/PheWeb/
```

Those files are not present in this cleaned repo. The Fig1 rebuild does not
require them. Reintroduce PheWeb only as a separate optional analysis if the
input files and coordinate conventions are restored.

