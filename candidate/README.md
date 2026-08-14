# Candidate sequence

`candidate_sequence.fasta` is the complete amino-acid sequence of the in silico multi-epitope candidate analyzed in the manuscript.

Release integrity is checked against `results/construct_summary.json`:

- expected length: 460 amino acids
- expected amino-acid sequence SHA-256: `c4f45f77ac6093a3fcea74d8a2af82a41a4b86720ecc31edb3870661dde8ba7a`

The validator computes the digest from the concatenated FASTA sequence, excluding the header and line breaks.
