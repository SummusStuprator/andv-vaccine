#!/usr/bin/env python3
"""Build a compact machine-readable summary of the released computational study."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
DATA = ROOT / "data"
CAND = ROOT / "candidate" / "candidate_sequence.fasta"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def csv_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fasta_sequence(path: Path) -> str:
    lines=[]
    for line in path.read_text(encoding="utf-8").splitlines():
        line=line.strip()
        if not line or line.startswith(">"):
            continue
        lines.append(line)
    return "".join(lines)


def main():
    hla = read_json(DATA / "hla_allele_panel.json")["hla_allele_panel"]
    raw = read_json(RES / "raw_predictions" / "deposit_manifest.json")
    cov = read_json(RES / "population_coverage.json")
    bcell = read_json(RES / "bcell_module.json")
    csum = read_json(RES / "construct_summary.json")
    c1 = csv_rows(RES / "cores_class1_scored.csv")
    c2 = csv_rows(RES / "cores_class2_scored.csv")
    s1 = csv_rows(RES / "final_epitopes_class1.csv")
    s2 = csv_rows(RES / "final_epitopes_class2.csv")
    sim1 = read_json(RES / "human_similarity_class1.json")
    sim2 = read_json(RES / "human_similarity_class2.json")

    if not CAND.is_file():
        raise FileNotFoundError(f"candidate sequence missing: {CAND.relative_to(ROOT)}")
    seq=fasta_sequence(CAND)
    seq_digest=hashlib.sha256(seq.encode("ascii")).hexdigest()

    def sim_counts(d):
        vals=[v for k,v in d.items() if k != "_meta"]
        out={"exact_match":0,"one_mismatch":0,"no_match":0}
        for v in vals:
            out[v.get("similarity_class","no_match")]+=1
        return out

    raw_c1 = sum(x["rows"] for x in raw["files"] if x["file"].startswith("class1_"))
    raw_c2 = sum(x["rows"] for x in raw["files"] if x["file"].startswith("class2_"))

    master = {
        "project": "ANDV computational vaccine analysis",
        "author": "Summus Stuprator",
        "scope": "computational analysis of an in silico vaccine candidate; no wet-lab experiments were performed",
        "data_provenance": {
            "genbank_M_segment": "MN258217.1",
            "genbank_S_segment": "MN258242.1",
            "canonical_sequence_file": "data/canonical_sequences.fasta",
            "native_reference_file": "data/native_genbank_reference.fasta",
            "data_dictionary": "data/data_dictionary.json",
        },
        "hla_panel": {
            "class_I_count": hla["class_I"]["count"],
            "class_II_count": hla["class_II"]["count"],
            "population_samples": 2693,
            "typed_loci": ["HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DQB1"],
        },
        "prediction_data": {
            "class_I_rows": raw_c1,
            "class_II_rows": raw_c2,
            "companion_archive": raw["companion_archive"],
            "companion_archive_sha256": raw["archive_sha256"],
            "dataset_doi": "10.5281/zenodo.21939568",
            "manifest": "results/raw_predictions/deposit_manifest.json",
        },
        "core_analysis": {
            "class_I_cores": len(c1),
            "class_II_cores": len(c2),
            "selected_class_I": len(s1),
            "selected_class_II": len(s2),
            "class_I_file": "results/cores_class1_scored.csv",
            "class_II_file": "results/cores_class2_scored.csv",
            "selected_class_I_file": "results/final_epitopes_class1.csv",
            "selected_class_II_file": "results/final_epitopes_class2.csv",
        },
        "human_proteome_similarity": {
            "class_I": sim_counts(sim1),
            "class_II": sim_counts(sim2),
            "terminology": "short-sequence similarity classes",
        },
        "population_coverage": cov,
        "bcell_prediction_summary": {
            "terminology_note": bcell.get("terminology_note"),
            "bepipred_linear_region_counts": {k: len(v) for k,v in bcell["bepipred_linear_regions"].items()},
            "discotope_high_scoring_residue_counts": {k: v["n_high_confidence_residues"] for k,v in bcell["discotope"].items()},
            "source": "archived external predictor tables; see docs/REPRODUCIBILITY.md",
        },
        "construct_summary": {
            **csum,
            "sequence_path": "candidate/candidate_sequence.fasta",
            "observed_sequence_sha256": seq_digest,
        },
        "release_content": "docs/RELEASE_CONTENT.md",
        "key_file_sha256": {
            "candidate/candidate_sequence.fasta": sha256(CAND),
            "data/canonical_sequences.fasta": sha256(DATA / "canonical_sequences.fasta"),
            "data/hla_allele_panel.json": sha256(DATA / "hla_allele_panel.json"),
            "results/cores_class1_scored.csv": sha256(RES / "cores_class1_scored.csv"),
            "results/cores_class2_scored.csv": sha256(RES / "cores_class2_scored.csv"),
            "results/population_coverage.json": sha256(RES / "population_coverage.json"),
        },
    }
    out = RES / "andv_vaccine_masterfile.json"
    out.write_text(json.dumps(master, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
