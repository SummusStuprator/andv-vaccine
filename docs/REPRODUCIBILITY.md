# Reproducing the analysis

The repository contains the inputs and downstream code used for the reported computational analyses, together with archived outputs from external services.

## Local calculations

The released code covers conservation summaries, overlap collapse and scoring, selected-set bookkeeping, genotype-based population coverage, human-proteome similarity post-processing, B-cell table post-processing, master-file generation, figure rendering, and release validation.

The large raw peptide-prediction CSVs are distributed through the linked Zenodo dataset (`10.5281/zenodo.21939568`) rather than stored directly in Git history. Their filenames, row counts, sizes, and hashes are recorded in `results/raw_predictions/deposit_manifest.json`.

## Candidate sequence

The complete 460-residue in silico candidate sequence is included at `candidate/candidate_sequence.fasta`. Its amino-acid sequence SHA-256 is recorded in `results/construct_summary.json`. The release validator checks both the sequence length and digest.

## External services

Prediction drivers for IEDB, NCBI sequence retrieval, and the 1000 Genomes population table are included. External services and databases can change, so a later rerun need not be byte-identical to the archived outputs used for the paper. Retrieval metadata are preserved with the corresponding files.

## B-cell predictor outputs

BepiPred and DiscoTope per-residue CSVs are archived in `results/`. The source structures originally used for the DiscoTope runs were not present in the audited working archive. The preserved tables support the downstream analysis reported in the paper, but not a fresh reconstruction of the structure-dependent prediction from primary structure files.

## HLA panel

The prediction panel contains 55 class-I and 42 class-II entries. Its metadata describe the frequency-based selection procedure used to form the panel. The panel-construction script itself was not retained, so `data/hla_allele_panel.json` is the released starting point for rerunning prediction stages.

## Population coverage

Coverage is calculated from the strong-binding HLA alleles of the selected representative peptides. Alleles contributed only by other members of the same overlap cluster are not credited.

The 1000 Genomes HLA table directly types HLA-A, -B, -C, -DRB1, and -DQB1. Class-II coverage is therefore reported twice: once using DRB1 only, and once using DRB1 plus the typed DQB1 beta chain of predicted DQ heterodimers.

## Release checks

`make release-check` rebuilds `MANIFEST.sha256` and runs the static release validator. These checks are intended to catch packaging and data-integrity mistakes; they are separate from biological interpretation of the candidate.
