# AR Cascade Perspective Figures

Reproducible code for figures in the AR cascade perspective manuscript. The original workflow below generates the allele-frequency figure (Figure 1); the Figure 4 carrier-evidence wheel and detailed supporting data are in [`figure4_carrier_evidence/`](figure4_carrier_evidence/).

The cleaned workflow focuses on the pieces that are present in this project now:
GWAS Catalog significant associations, ClinVar variants, and gnomAD v4 allele
frequencies. Older PheWeb-dependent analysis is documented as legacy context but
is not required to rebuild Fig1 from this repo.

## Repository Layout

```text
.
|-- scripts/
|   |-- extract_gnomad_af.sh
|   |-- prep_gnomad_clinvar_join.py
|   `-- make_fig1_gwas_clinvar.R
|-- data/
|   |-- raw/
|   |-- interim/
|   `-- processed/
|-- results/
|   `-- figures/
|-- docs/
|   |-- data_sources.md
|   `-- source_script_manifest.md
|-- environment.yml
|-- Makefile
`-- README.md
```

## Inputs

Large data files are intentionally not tracked in git. Place or symlink inputs
into these locations:

```text
data/raw/gwas_catalog_v1.0.2-associations_e115_r2026-04-27_full.tsv
data/raw/variant_summary.txt.gz
data/interim/gnomad_v4_AF/gnomad_v4_AF_chr1.txt.gz
...
data/interim/gnomad_v4_AF/gnomad_v4_AF_chrY.txt.gz
```

The current working folder used lab-local source paths for GWAS Catalog release
`e115 r2026-04-27`, ClinVar `May 2026`, and gnomAD v4. See
`docs/data_sources.md` for details.

## Quick Start

Create the environment, or install equivalent tools manually:

```bash
conda env create -f environment.yml
conda activate ar-cascade-fig1
```

If the per-chromosome gnomAD AF files already exist, run:

```bash
make prep
make figure
```

To stream per-chromosome gnomAD AF files from UCSC first:

```bash
make extract-gnomad
```

On SLURM, the extractor can also run as an array job:

```bash
sbatch --array=1-24 scripts/extract_gnomad_af.sh
```

## Main Commands

Prepare ClinVar and GWAS Catalog joins against gnomAD:

```bash
python scripts/prep_gnomad_clinvar_join.py \
  --gnomad-dir data/interim/gnomad_v4_AF \
  --clinvar data/raw/variant_summary.txt.gz \
  --gwas-catalog data/raw/gwas_catalog_v1.0.2-associations_e115_r2026-04-27_full.tsv \
  --outdir data/processed
```

Generate Fig1:

```bash
Rscript scripts/make_fig1_gwas_clinvar.R \
  --gwas-catalog data/raw/gwas_catalog_v1.0.2-associations_e115_r2026-04-27_full.tsv \
  --clinvar data/raw/variant_summary.txt.gz \
  --joined-clinvar-gnomad data/processed/clinvar_gnomad_v4_joined_allele.tsv \
  --outdir results/figures
```

Outputs:

```text
results/figures/fig1_gwas_catalog_clinvar_gnomad.pdf
results/figures/fig1_gwas_catalog_clinvar_gnomad.png
results/figures/fig1_gwas_catalog_clinvar_gnomad_summary.tsv
```

## Notes For The Manuscript

The original v3 plotting script computed both `AF` and `AF_GRPMAX` from the
gnomAD join, but plotted `AF`. This cleaned script keeps `AF` as the default to
match the current figure-producing code. If the figure legend should use
grpmax, rerun with:

```bash
Rscript scripts/make_fig1_gwas_clinvar.R --clinvar-af-column AF_GRPMAX
```

The PheWeb overlap workbook from the older script is not reproduced here because
the required PheWeb result files are not part of the current cleaned workflow.

Before public release, choose and add the appropriate license for the manuscript
code.
