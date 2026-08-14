# Provenance

## Canonical ANDV inputs

The analysis uses GenBank accessions `MN258217.1` for the M segment/GPC reference and `MN258242.1` for the S segment/nucleocapsid reference. Checksums and stored-sequence descriptions are recorded in `data/data_dictionary.json` and `data/canonical_sequences.fasta`.

## HLA population data

The genotype model uses the 2,693-sample 1000 Genomes HLA typing panel described by Abi-Rached et al. (2018). The source table and its upstream manifest are under `data/population/`. The genotype table types HLA-A, -B, -C, -DRB1, and -DQB1; coverage outside those typed loci is not inferred by the local genotype model.

## Conservation data

`data/conservation/` contains the sequence panels, fetch manifest, accession information, and checksums used for the conservation analysis. MAFFT v7.526 is the recorded alignment dependency.

## Prediction data

The Git repository keeps prediction job manifests and per-chunk metadata. The large raw prediction CSVs are packaged separately. `results/raw_predictions/deposit_manifest.json` records their row counts and SHA-256 hashes and the SHA-256 of the companion archive.

## External predictor outputs

BepiPred and DiscoTope per-residue outputs are preserved in `results/`. Their downstream processing is reproducible from those tables. The upstream structure files used for the DiscoTope run were not available in the audited workspace and are not represented as if they were.

## Candidate sequence

The released candidate is stored at `candidate/candidate_sequence.fasta`. `results/construct_summary.json` records the expected 460-residue length and amino-acid sequence SHA-256 used by the release validator.
