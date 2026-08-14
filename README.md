# ANDVOR - Hantavirus Andes Open Research Vaccine

This repository contains the computational research package for an in silico multi-epitope Andes virus (ANDV; *Orthohantavirus andesense*) vaccine candidate designated "ANDVOR".

**Author:** Summus Stuprator

The project combines peptide-HLA prediction, reduction of overlapping predictions to non-redundant regions, sequence-conservation analysis, human-proteome sequence comparison, genotype-based population-coverage modeling, secondary B-cell prediction, and candidate-level sequence descriptors. This is strictly an in-silico vaccine candidate, meaning that all designs were done only on computer systems. None of the content here has been validated in real-world laboratory conditions at this time. 
## Main results

We analyzed 342,430 class-I and 64,680 class-II peptide-allele prediction rows and succesfully reduced strong-binding predictions to 836 unique class-I peptides and 308 unique class-II peptides.

Overlap collapse produced 221 class-I and 73 class-II non-redundant cores. The candidate set contains 15 class-I and 10 class-II representatives.

Genotype-model combined coverage was 97.55% under the DRB1+DQB1 typed-loci model and 93.50% in the DRB1-only sensitivity model. Conservation, human-proteome sequence similarity, B-cell prediction outputs, and aggregate candidate properties are included with the release.

## Repository contents

```text
candidate/             Complete in silico candidate sequence and candidate metadata
data/                  Viral sequence, HLA, conservation, and population inputs
pipeline/              Analysis code
results/               Derived results and archived predictor outputs
results/raw_predictions/
                       Metadata and manifest for the large raw prediction archive
paper/                 Manuscript, HTML version, figures, and source data
paper/scripts/         Figure-generation code
docs/                  Data access, provenance, and analysis notes
.github/               Validation workflow and contribution templates
```

The complete candidate amino-acid sequence is stored at `candidate/candidate_sequence.fasta`. `results/construct_summary.json` records its expected length, SHA-256 sequence digest, and aggregate sequence-derived properties.

## Raw prediction dataset

The large raw prediction CSVs are archived separately on Zenodo rather than committed to Git history.

**Dataset DOI:** `10.5281/zenodo.21939568`

`results/raw_predictions/deposit_manifest.json` records the archive filename, per-file row counts and hashes, and the SHA-256 of the companion archive.

## Manuscript

- `paper/ANDV_manuscript_v1.0.pdf` - publication PDF
- `paper/ANDV_manuscript_v1.0.docx` - editable manuscript
- `paper/manuscript.md` - Markdown source
- `paper/index.html` - browser-readable HTML version

The HTML file can also be served with GitHub Pages.

## Validation

Install the Python requirements and run:

```bash
python -m pip install -r requirements.txt
make release-check
```

The release validator checks authorship metadata, candidate-sequence integrity, expected prediction-row totals, corrected population-coverage outputs, stale release-language markers, and `MANIFEST.sha256`.

To regenerate the six manuscript figures and their machine-readable source-data CSVs:

```bash
make figures
```

Some upstream analyses depend on external services. Retrieval metadata and archived outputs identify the data state used for the paper.

## Citation and licensing

Citation metadata are in `CITATION.cff`. Code is licensed under Apache-2.0. Original project documentation, figures, and project-created research data are offered under CC BY 4.0 unless a file states otherwise. Imported datasets and third-party material retain their original terms; see `THIRD_PARTY.md` and `NOTICE.md`.
