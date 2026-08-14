#!/usr/bin/env python3
"""Static checks for the versioned research release."""
from __future__ import annotations
import hashlib
import json
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ERRORS=[]
TEXT_SUFFIXES={'.md','.txt','.json','.py','.sh','.cff','.yml','.yaml','.csv','.html'}

patterns={
    'absolute Windows user path': re.compile(r'C:\\\\Users\\\\', re.I),
    'placeholder owner': re.compile(r'<owner>', re.I),
    'placeholder ORCID': re.compile(r'0000-0000-0000-0000'),
    'internal version-history label': re.compile(r'\bv2\.[01](?:\.0)?\b', re.I),
    'legacy/review-response wording': re.compile(r'\b(?:legacy construct|response_to_review|superseded by|previous version)\b', re.I),
    'stale construct-exclusion wording': re.compile(
        r'(?:exact engineered construct sequence is (?:outside|not part)|'
        r'public Git tree does not contain the exact engineered construct|'
        r'exact construct material outside the public|'
        r'keep exact construct material outside)', re.I),
}
for p in ROOT.rglob('*'):
    if not p.is_file() or p == Path(__file__).resolve() or '.git' in p.parts or '__pycache__' in p.parts or p.suffix.lower() not in TEXT_SUFFIXES:
        continue
    try:
        text=p.read_text(encoding='utf-8')
    except UnicodeDecodeError:
        continue
    for label,pat in patterns.items():
        if pat.search(text):
            ERRORS.append(f"{label}: {p.relative_to(ROOT)}")

# Authorship metadata.
if 'Summus Stuprator' not in (ROOT/'AUTHORS.md').read_text(encoding='utf-8'):
    ERRORS.append('AUTHORS.md does not identify Summus Stuprator')
cff=(ROOT/'CITATION.cff').read_text(encoding='utf-8')
if 'family-names: Stuprator' not in cff or 'given-names: Summus' not in cff:
    ERRORS.append('CITATION.cff author metadata mismatch')
if 'version: 1.0.0' not in cff:
    ERRORS.append('CITATION.cff release version is not 1.0.0')

# Candidate sequence integrity.
csum=json.loads((ROOT/'results/construct_summary.json').read_text(encoding='utf-8'))
candidate=ROOT/csum.get('sequence_path','candidate/candidate_sequence.fasta')
if not candidate.is_file():
    ERRORS.append(f'candidate sequence missing: {candidate.relative_to(ROOT)}')
else:
    lines=[]
    headers=0
    for raw in candidate.read_text(encoding='utf-8').splitlines():
        line=raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            headers += 1
            continue
        lines.append(line)
    seq=''.join(lines).upper()
    if headers < 1:
        ERRORS.append('candidate FASTA has no header')
    if not re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWY]+', seq):
        ERRORS.append('candidate FASTA contains non-standard amino-acid characters')
    if len(seq) != int(csum.get('length_aa',-1)):
        ERRORS.append(f"candidate length {len(seq)} != {csum.get('length_aa')}")
    got=hashlib.sha256(seq.encode('ascii')).hexdigest()
    if got != csum.get('sha256'):
        ERRORS.append(f"candidate sequence SHA-256 {got} != construct_summary sha256")
if csum.get('sequence_publicly_included') is not True:
    ERRORS.append('construct_summary does not mark candidate sequence as included')

# Raw prediction archive manifest.
dep=json.loads((ROOT/'results/raw_predictions/deposit_manifest.json').read_text(encoding='utf-8'))
c1=sum(x['rows'] for x in dep['files'] if x['file'].startswith('class1_'))
c2=sum(x['rows'] for x in dep['files'] if x['file'].startswith('class2_'))
if c1 != 342430:
    ERRORS.append(f'class-I raw row count {c1} != 342430')
if c2 != 64680:
    ERRORS.append(f'class-II raw row count {c2} != 64680')
if not re.fullmatch(r'[0-9a-f]{64}', dep.get('archive_sha256','')):
    ERRORS.append('invalid raw archive SHA-256')

# Corrected population-coverage results.
cov=json.loads((ROOT/'results/population_coverage.json').read_text(encoding='utf-8'))
checks=[
    ('typed-loci combined', cov['drb1_dqb1']['combined_pct'], 97.55),
    ('typed-loci class II', cov['drb1_dqb1']['class_II_pct'], 97.99),
    ('DRB1-only combined', cov['drb1_only']['combined_pct'], 93.50),
    ('DRB1-only class II', cov['drb1_only']['class_II_pct'], 93.95),
]
for label,got,want in checks:
    if abs(float(got)-want) > 1e-9:
        ERRORS.append(f'{label} {got} != {want}')

# Canonical result names only.
for p in (ROOT/'results').glob('*v21*'):
    ERRORS.append(f"version-suffixed result remains: {p.name}")
for p in (ROOT/'results').glob('*v2*'):
    ERRORS.append(f"version-suffixed result remains: {p.name}")

# Figure and source-data inventory.
for stem in [
    'fig1_workflow',
    'fig2_analysis_reduction',
    'fig3_selected_conservation',
    'fig4_representative_allele_breadth',
    'fig5_coverage_by_superpopulation',
    'fig6_candidate_architecture',
]:
    if not (ROOT/'paper/figures'/f'{stem}.png').is_file():
        ERRORS.append(f'missing figure: paper/figures/{stem}.png')
    if not (ROOT/'paper/source_data'/f'{stem}.csv').is_file():
        ERRORS.append(f'missing figure source data: paper/source_data/{stem}.csv')

# Manifest integrity.
mf=ROOT/'MANIFEST.sha256'
if not mf.exists():
    ERRORS.append('MANIFEST.sha256 missing')
else:
    for n,line in enumerate(mf.read_text(encoding='utf-8').splitlines(),1):
        if not line.strip():
            continue
        try:
            digest,rel=line.split('  ',1)
        except ValueError:
            ERRORS.append(f'bad manifest line {n}')
            continue
        p=ROOT/rel
        if not p.is_file():
            ERRORS.append(f'manifest path missing: {rel}')
            continue
        got=hashlib.sha256(p.read_bytes()).hexdigest()
        if got != digest:
            ERRORS.append(f'manifest mismatch: {rel}')

if ERRORS:
    for e in ERRORS:
        print('FAIL',e)
    raise SystemExit(1)
print('PASS release integrity checks')
