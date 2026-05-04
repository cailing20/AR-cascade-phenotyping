GWAS_CATALOG ?= data/raw/gwas_catalog_v1.0.2-associations_e115_r2026-04-27_full.tsv
CLINVAR ?= data/raw/variant_summary.txt.gz
GNOMAD_AF_DIR ?= data/interim/gnomad_v4_AF
PROCESSED_DIR ?= data/processed
FIGURE_DIR ?= results/figures
FIGURE_PREFIX ?= fig1_gwas_catalog_clinvar_gnomad

.PHONY: help extract-gnomad prep figure all clean

help:
	@echo "Targets:"
	@echo "  make extract-gnomad  Stream gnomAD v4 AF fields into per-chromosome files"
	@echo "  make prep            Join ClinVar and GWAS Catalog records to gnomAD AF files"
	@echo "  make figure          Build the GWAS Catalog vs ClinVar/gnomAD figure"
	@echo "  make all             Run prep and figure"
	@echo ""
	@echo "Override paths like:"
	@echo "  make figure GWAS_CATALOG=/path/gwas.tsv CLINVAR=/path/variant_summary.txt.gz"

extract-gnomad:
	OUTDIR="$(GNOMAD_AF_DIR)" bash scripts/extract_gnomad_af.sh

prep:
	python scripts/prep_gnomad_clinvar_join.py \
		--gnomad-dir "$(GNOMAD_AF_DIR)" \
		--clinvar "$(CLINVAR)" \
		--gwas-catalog "$(GWAS_CATALOG)" \
		--outdir "$(PROCESSED_DIR)"

figure:
	Rscript scripts/make_fig1_gwas_clinvar.R \
		--gwas-catalog "$(GWAS_CATALOG)" \
		--clinvar "$(CLINVAR)" \
		--joined-clinvar-gnomad "$(PROCESSED_DIR)/clinvar_gnomad_v4_joined_allele.tsv" \
		--outdir "$(FIGURE_DIR)" \
		--prefix "$(FIGURE_PREFIX)"

all: prep figure

clean:
	@echo "Generated data are ignored by git. Remove data/interim, data/processed, or results manually if needed."

