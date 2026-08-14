# Data access

| Data class | Location | Format | Notes |
|---|---|---|---|
| Complete candidate sequence | `candidate/candidate_sequence.fasta` | FASTA | 460-residue in silico candidate; digest in `results/construct_summary.json` |
| Candidate summary | `results/construct_summary.json` | JSON | Sequence-derived descriptors and candidate identifier |
| Canonical viral sequences | `data/canonical_sequences.fasta` | FASTA | Stored analysis sequences |
| Native GenBank references | `data/native_genbank_reference.fasta` | FASTA | Reference translations used for coordinate mapping |
| HLA prediction panel | `data/hla_allele_panel.json` | JSON | 55 class-I and 42 class-II entries |
| 1000 Genomes HLA types | `data/population/` | TSV/text | 2,693 samples |
| Conservation panels | `data/conservation/` | FASTA/JSON | Sequence panels and accession manifest |
| Raw peptide-HLA predictions | Zenodo DOI `10.5281/zenodo.21939568` | ZIP / CSV / JSON | 342,430 class-I and 64,680 class-II prediction rows |
| Raw archive manifest | `results/raw_predictions/deposit_manifest.json` | JSON | Per-file rows, sizes, hashes, and archive digest |
| Core and selected-set tables | `results/` | CSV | Non-redundant cores and selected representatives |
| Human-proteome comparison | `results/human_similarity_class*.json` | JSON | Sequence-similarity results |
| Genotype coverage | `results/population_coverage.json` | JSON | Overall and superpopulation estimates with bootstrap intervals |
| B-cell predictor outputs | `results/bepipred_*.csv`, `results/discotope_*.csv` | CSV | Archived external predictor outputs |
| Manuscript | `paper/ANDV_manuscript_v1.0.pdf`, `.docx`, `paper/manuscript.md`, `paper/index.html` | PDF/DOCX/Markdown/HTML | Publication and browser-readable versions |
| Figure source data | `paper/source_data/` | CSV | Machine-readable values for all six manuscript figures |
| Release inventory | `MANIFEST.sha256` | text | SHA-256 over tracked release files |
