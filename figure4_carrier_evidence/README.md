# Carrier evidence wheel

Generate the current Section 7, Figure 4 carrier-phenotype evidence wheel from the bundled detailed evidence workbook. This is a **preliminary research figure**: source review, inheritance adjudication, and publication-quality review are still needed. It is not a systematic-review result, and blank cells are not evidence of absent phenotypes.

## Files

- `Supplementary_Table_Figure_4_carrier_evidence.xlsx` is the detailed evidence and plotting workbook bundled in this folder. The `Figure marks` sheet is the **direct plotting input**. Each row describes one comparison mark, including the original mark ID and evidence IDs that link back to `Evidence details` and `Papers`.
- `generate_wheel.py` reads the workbook and creates a vector PDF and SVG, a high-resolution PNG, and a count summary.
- `config/figure.json` defines domain order, gene order, ring order, colors, and size-bin geometry.

The original sheets and their held-for-review records remain in the workbook. The `Figure marks` sheet captures the current 100-mark snapshot, including source and allele-context flags. When changing an observation in `Evidence details`, **also review and update its corresponding row or rows in `Figure marks`**: changes to `Evidence details`, `Papers`, or `Held for review` are not propagated into the plot automatically. The plotting script does not screen or adjudicate sources. When changing the number of marks, publications, or genes, update the Figure 4 legend and manuscript as appropriate.

## Setup

Requires Python 3.10 or newer. In a terminal, from this folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, activate the virtual environment with `.venv\\Scripts\\activate` instead.

## Generate

```bash
python generate_wheel.py
```

By default, the script reads the bundled workbook and writes to `outputs/`:

- `carrier_study_wheel.pdf` — vector artwork
- `carrier_study_wheel.svg` — vector artwork, including accessible mark labels
- `carrier_study_wheel.png` — high-resolution image
- `summary.json` — counts and the workbook path

For a different input workbook or output folder:

```bash
python generate_wheel.py --workbook /path/to/table.xlsx --output-dir /path/to/output
```

Use `python generate_wheel.py --help` for options. Generated outputs are ignored by Git so reruns do not clutter version control. The workbook is versioned in this folder. To use a local working copy, specify `--workbook`.

## Editing the plotting data

In `Figure marks`, `direction` must be `Positive`, `Null`, `Protective`, or `Mixed`; `carrier_n_max` is the **endpoint-specific carrier/donor N actually used for size**, and must be a positive whole number or blank if unknown. Do not enter a cohort enrollment count as an endpoint count unless the source supports it. Marker sizes are computed from that column on each run. `carrier_n_min` is optional contextual information; if retained, it must not exceed `carrier_n_max`. When correcting `carrier_n_max`, also correct or clear `carrier_n_min` if needed. `abstract_only` and `qualified_context` accept `Yes` or `No`; the latter supplies the dagger and is not a complete inheritance annotation. Use the linked evidence records and original sources to review allele, phenotype, cohort, and overlap before changes.

To include a newly curated comparison, add a uniquely identified row with its gene, domain, publication, measurement category, source flags, and evidence links. A gene must be placed in its intended domain in `config/figure.json`; the script reports an error for unknown genes or domain mismatches. Add/remove rows only after reviewing the source and update the `Papers` and `Evidence details` sheets and the figure legend accordingly. No script can infer whether a finding should be included from gene name alone.

The five rings represent measurement categories rather than effect size, severity, dominance, or clinical actionability. Multiple marks or papers may share participants and are not necessarily independent replications. Some included evidence is abstract-only, and exact variant notation for the *IL2RB* source remains unresolved in the workbook.

## Reproduction check

The supplied baseline workbook produces 33 genes, 70 publications, and 100 marks. On the bundled environment used to prepare this folder, its SVG and PDF were byte-for-byte identical to the prior Figure 4 artwork. PNGs can vary slightly by rendering engine. A test on a disposable workbook changed one mark from a 173-person positive finding to a 7-person null finding; rerunning produced the corresponding smaller, unfilled symbol. The bundled workbook was not modified by that test.
