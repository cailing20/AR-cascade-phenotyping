# Source Script Manifest

This folder is a cleaned GitHub-ready extraction of the manuscript working
folder at:

```text
AR_cascade_perspective/
```

## Cleaned Scripts

| Clean script | Source script(s) | What changed |
| --- | --- | --- |
| `scripts/extract_gnomad_af.sh` | `scripts/gnomad_prep.sh`, `scripts/gnomad_extract.sh`, `scripts/gnomad_extract_retry.sh` | Consolidated local and SLURM-array extraction; removed hard-coded output paths; kept `bcftools query` logic for `AF` and `AF_grpmax`. |
| `scripts/prep_gnomad_clinvar_join.py` | `scripts/prep_gnomad_clinvar_join_allele_v2.py` | Added command-line arguments, path defaults inside the repo, stricter input-column checks, and CSV writer output. |
| `scripts/make_fig1_gwas_clinvar.R` | `scripts/GWAS_significant_association_count_by_MAF_v3.R` | Removed hard-coded lab paths, removed required PheWeb dependency, made the plotted ClinVar AF column configurable, and writes PDF, PNG, and a small summary TSV. |

## Legacy Material Not Carried Forward

The source folder also contains:

```text
scripts/prior/assemble_IEM_gene_phenotypes.R
scripts/prior/GWAS_significant_association_count_by_MAF.R
```

Those appear to belong to earlier phenotype/PheWeb iterations and are not
needed for the current GWAS Catalog plus ClinVar/gnomAD Fig1 workflow. Keep them
in the manuscript working folder for provenance, but do not upload them as part
of the cleaned repo unless a separate legacy-analysis branch is needed.

