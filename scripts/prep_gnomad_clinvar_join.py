#!/usr/bin/env python3
"""Prepare ClinVar and GWAS Catalog joins against per-chromosome gnomAD AF files.

The gnomAD input files are expected to be tab-separated, gzip-compressed files
with columns:

    CHROM POS REF ALT AF AF_GRPMAX

They can be generated with scripts/extract_gnomad_af.sh.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, Iterator, List, Tuple


DEFAULT_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
ClinVarMeta = Dict[str, str]


def parse_chromosomes(value: str) -> List[str]:
    if value.lower() == "all":
        return DEFAULT_CHROMS

    chroms: List[str] = []
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            start, end = token.split("-", 1)
            chroms.extend(str(i) for i in range(int(start), int(end) + 1))
        else:
            chroms.append(token.replace("chr", ""))
    return chroms


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("rt", newline="")


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.startswith("chr"):
        return chrom
    return f"chr{chrom}"


def require_columns(reader: csv.DictReader, columns: Iterable[str], path: Path) -> None:
    missing = [col for col in columns if col not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")


def load_clinvar(
    clinvar_path: Path,
) -> Tuple[
    DefaultDict[Tuple[str, str, str, str], List[ClinVarMeta]],
    DefaultDict[Tuple[str, str], List[ClinVarMeta]],
]:
    """Load ClinVar GRCh38 records into allele-aware and position-only maps."""

    allele_dict: DefaultDict[Tuple[str, str, str, str], List[ClinVarMeta]]
    pos_dict: DefaultDict[Tuple[str, str], List[ClinVarMeta]]
    allele_dict = defaultdict(list)
    pos_dict = defaultdict(list)

    with open_text(clinvar_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(
            reader,
            [
                "Assembly",
                "Chromosome",
                "Start",
                "ReferenceAllele",
                "AlternateAllele",
                "ClinicalSignificance",
                "GeneSymbol",
                "Name",
                "PhenotypeList",
            ],
            clinvar_path,
        )

        for row in reader:
            if row["Assembly"] != "GRCh38":
                continue

            chrom = normalize_chrom(row["Chromosome"])
            pos = row["Start"].strip()
            ref = row["ReferenceAllele"].strip().upper()
            alt = row["AlternateAllele"].strip().upper()
            meta = {
                "ClinicalSignificance": row["ClinicalSignificance"],
                "GeneSymbol": row["GeneSymbol"],
                "Name": row["Name"],
                "PhenotypeList": row["PhenotypeList"],
            }

            if not pos:
                continue

            if ref in {"", "NA"} or alt in {"", "NA"}:
                pos_dict[(chrom, pos)].append(meta)
                continue

            for alt_allele in alt.split(","):
                alt_allele = alt_allele.strip()
                if alt_allele and alt_allele != "NA":
                    allele_dict[(chrom, pos, ref, alt_allele)].append(meta)
                else:
                    pos_dict[(chrom, pos)].append(meta)

    return allele_dict, pos_dict


def iter_gnomad_records(path: Path) -> Iterator[Tuple[str, str, str, List[str], List[str], List[str]]]:
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue

            chrom, pos, ref, alt_string, af_string, grpmax_string = parts[:6]
            alts = [alt.strip().upper() for alt in alt_string.split(",")]
            afs = af_string.split(",")
            grpmaxs = grpmax_string.split(",")

            while len(afs) < len(alts):
                afs.append(".")
            while len(grpmaxs) < len(alts):
                grpmaxs.append(".")

            yield chrom, pos, ref.upper(), alts, afs, grpmaxs


def write_clinvar_gnomad_join(
    gnomad_dir: Path,
    outdir: Path,
    chromosomes: List[str],
    allele_dict: DefaultDict[Tuple[str, str, str, str], List[ClinVarMeta]],
    pos_dict: DefaultDict[Tuple[str, str], List[ClinVarMeta]],
) -> Path:
    out_path = outdir / "clinvar_gnomad_v4_joined_allele.tsv"
    header = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "AF",
        "AF_GRPMAX",
        "ClinicalSignificance",
        "GeneSymbol",
        "Name",
        "PhenotypeList",
        "match_type",
    ]

    matched_allele = 0
    matched_pos = 0

    with out_path.open("w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)

        for chrom_num in chromosomes:
            chrom = f"chr{chrom_num}"
            gnomad_path = gnomad_dir / f"gnomad_v4_AF_chr{chrom_num}.txt.gz"
            if not gnomad_path.exists():
                print(f"WARNING: {gnomad_path} not found; skipping {chrom}")
                continue

            print(f"Processing ClinVar x gnomAD: {chrom}")
            for g_chrom, g_pos, g_ref, alts, afs, grpmaxs in iter_gnomad_records(gnomad_path):
                for idx, (alt, af, grpmax) in enumerate(zip(alts, afs, grpmaxs)):
                    allele_key = (g_chrom, g_pos, g_ref, alt)
                    for cv in allele_dict.get(allele_key, []):
                        writer.writerow(
                            [
                                g_chrom,
                                g_pos,
                                g_ref,
                                alt,
                                af,
                                grpmax,
                                cv["ClinicalSignificance"],
                                cv["GeneSymbol"],
                                cv["Name"],
                                cv["PhenotypeList"],
                                "allele",
                            ]
                        )
                        matched_allele += 1

                    # Position-only ClinVar records do not identify the ALT allele.
                    # Write them once per gnomAD site to avoid multiplying rows.
                    if idx == 0:
                        pos_key = (g_chrom, g_pos)
                        for cv in pos_dict.get(pos_key, []):
                            writer.writerow(
                                [
                                    g_chrom,
                                    g_pos,
                                    g_ref,
                                    ",".join(alts),
                                    af,
                                    grpmax,
                                    cv["ClinicalSignificance"],
                                    cv["GeneSymbol"],
                                    cv["Name"],
                                    cv["PhenotypeList"],
                                    "position_only",
                                ]
                            )
                            matched_pos += 1

    print("ClinVar x gnomAD join complete")
    print(f"  Allele-matched rows: {matched_allele:,}")
    print(f"  Position-only rows:  {matched_pos:,}")
    print(f"  Output: {out_path}")
    return out_path


def load_gwas_positions(
    gwas_catalog: Path,
    pvalue_mlog_threshold: float,
) -> DefaultDict[Tuple[str, str], List[Dict[str, str]]]:
    gwas_pos: DefaultDict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)

    with gwas_catalog.open("rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(
            reader,
            ["CHR_ID", "CHR_POS", "PVALUE_MLOG", "RISK ALLELE FREQUENCY", "DISEASE/TRAIT"],
            gwas_catalog,
        )

        for row in reader:
            try:
                pvalue_mlog = float(row["PVALUE_MLOG"])
            except ValueError:
                continue

            if pvalue_mlog <= pvalue_mlog_threshold:
                continue

            chrom = row["CHR_ID"].strip()
            pos = row["CHR_POS"].strip()
            if not chrom or ";" in chrom or not pos:
                continue

            gwas_pos[(normalize_chrom(chrom), pos)].append(
                {
                    "RISK_AF": row.get("RISK ALLELE FREQUENCY", ""),
                    "TRAIT": row.get("DISEASE/TRAIT", ""),
                    "PVAL_MLOG": row["PVALUE_MLOG"],
                }
            )

    return gwas_pos


def write_gwas_gnomad_join(
    gnomad_dir: Path,
    outdir: Path,
    chromosomes: List[str],
    gwas_pos: DefaultDict[Tuple[str, str], List[Dict[str, str]]],
) -> Path:
    out_path = outdir / "gwas_catalog_gnomad_v4_joined.tsv"
    header = ["CHROM", "POS", "REF", "ALT", "AF", "AF_GRPMAX", "RISK_AF", "TRAIT", "PVAL_MLOG"]
    matched = 0

    with out_path.open("w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)

        for chrom_num in chromosomes:
            chrom = f"chr{chrom_num}"
            gnomad_path = gnomad_dir / f"gnomad_v4_AF_chr{chrom_num}.txt.gz"
            if not gnomad_path.exists():
                continue

            print(f"Processing GWAS Catalog x gnomAD: {chrom}")
            for g_chrom, g_pos, g_ref, alts, afs, grpmaxs in iter_gnomad_records(gnomad_path):
                records = gwas_pos.get((g_chrom, g_pos))
                if not records:
                    continue

                first_af = afs[0] if afs else "."
                first_grpmax = grpmaxs[0] if grpmaxs else "."
                for gw in records:
                    writer.writerow(
                        [
                            g_chrom,
                            g_pos,
                            g_ref,
                            ",".join(alts),
                            first_af,
                            first_grpmax,
                            gw["RISK_AF"],
                            gw["TRAIT"],
                            gw["PVAL_MLOG"],
                        ]
                    )
                    matched += 1

    print("GWAS Catalog x gnomAD join complete")
    print(f"  Rows:   {matched:,}")
    print(f"  Output: {out_path}")
    return out_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gnomad-dir", type=Path, default=Path("data/interim/gnomad_v4_AF"))
    parser.add_argument("--clinvar", type=Path, default=Path("data/raw/variant_summary.txt.gz"))
    parser.add_argument(
        "--gwas-catalog",
        type=Path,
        default=Path("data/raw/gwas_catalog_v1.0.2-associations_e115_r2026-04-27_full.tsv"),
    )
    parser.add_argument("--outdir", type=Path, default=Path("data/processed"))
    parser.add_argument("--chromosomes", default="all", help="Comma list such as 1-22,X,Y or 'all'.")
    parser.add_argument("--pvalue-mlog-threshold", type=float, default=8.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    chromosomes = parse_chromosomes(args.chromosomes)

    print("Loading ClinVar GRCh38 records")
    allele_dict, pos_dict = load_clinvar(args.clinvar)
    print(f"  Allele-matchable keys: {len(allele_dict):,}")
    print(f"  Position-only keys:    {len(pos_dict):,}")

    write_clinvar_gnomad_join(
        args.gnomad_dir,
        args.outdir,
        chromosomes,
        allele_dict,
        pos_dict,
    )

    print("Loading GWAS Catalog significant positions")
    gwas_pos = load_gwas_positions(args.gwas_catalog, args.pvalue_mlog_threshold)
    print(f"  Positions: {len(gwas_pos):,}")

    write_gwas_gnomad_join(args.gnomad_dir, args.outdir, chromosomes, gwas_pos)


if __name__ == "__main__":
    main()

