# Contributing

Contributions are welcome for documentation, data corrections, analysis code, tests, figures, and reproducibility work.

## Ground rules

1. Quantitative changes should identify the input data and code path that produced them.
2. Preserve accession numbers, retrieval dates, hashes, and external-tool metadata when updating source material.
3. Keep computational predictions distinct from experimental measurements in documentation and discussion.
4. Do not silently replace files that belong to a published release. Use a new version or a clearly described change.
5. Run `make release-check` before opening a pull request.

## Pull requests

Describe what changed, how you checked it, and which result files are affected. If an external service was rerun, record the retrieval date and preserve the returned metadata.

## Reproducibility reports

For a mismatch, include the command used, expected result, observed result, relevant file or commit, and local environment. Do not post private identity information in public issues.
