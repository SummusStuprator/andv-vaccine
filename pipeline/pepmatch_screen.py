#!/usr/bin/env python3
"""Screen unique strong-binding peptides against the human proteome with PEPMatch.

The output labels exact_match, one_mismatch, and no_match are sequence-similarity classes used by this computational pipeline. They are not a validated clinical risk taxonomy.
"""
import glob
import json
import os
import time
from datetime import datetime, timezone

import pandas as pd
import requests

BASE = "https://api-nextgen-tools.iedb.org/api/v1"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
PROT = {"ANDV_Gn": 1, "ANDV_Gc": 2, "ANDV_N": 3}


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def unique_strong_peptides(cls):
    if cls == 1:
        df = pd.concat([pd.read_csv(f) for f in
                        sorted(glob.glob(os.path.join(RES, "raw_predictions/class1_chunk*.csv")))])
        s = df[df["netmhcpan_el_percentile"] < 0.5]
    else:
        df = pd.concat([pd.read_csv(f) for f in
                        sorted(glob.glob(os.path.join(RES, "raw_predictions/class2_chunk[0-9].csv")))])
        s = df[df["netmhciipan_el_percentile"] < 2.0]
    return sorted(s["peptide"].unique())


def submit(peptides):
    payload = {
        "pipeline_title": "andv_v21_pepmatch",
        "run_stage_range": [1, 1],
        "stages": [{
            "stage_number": 1,
            "tool_group": "pepmatch",
            "input_sequence_text": "\n".join(peptides),
            "input_parameters": {"proteome": "Human", "mismatch": 1, "best_match": True},
        }],
    }
    r = requests.post(f"{BASE}/pipeline", json=payload, timeout=120)
    r.raise_for_status()
    return r.json()


def poll(pipeline_id, interval=20, timeout_min=60):
    """Poll the stage endpoint until done; return the stage_result_id."""
    deadline = time.time() + timeout_min * 60
    while time.time() < deadline:
        r = requests.get(f"{BASE}/pipeline/{pipeline_id}", timeout=60)
        r.raise_for_status()
        stages = r.json().get("stages", [])
        if not stages:
            time.sleep(interval)
            continue
        surl = stages[0].get("stage_url")
        sd = requests.get(surl, timeout=60).json() if surl else stages[0]
        st = sd.get("stage_status")
        if st == "done":
            return sd.get("stage_result_id")
        if st in ("error", "failed"):
            raise RuntimeError(f"pipeline failed: {sd.get('stage_messages')}")
        time.sleep(interval)
    raise TimeoutError("pepmatch pipeline timed out")


def fetch_table(result_id):
    r = requests.get(f"{BASE}/results/{result_id}", timeout=300)
    r.raise_for_status()
    d = r.json()
    for b in d["data"]["results"]:
        if b.get("type") == "peptide_table":
            cols = [c["name"] for c in b["table_columns"]]
            return [dict(zip(cols, row)) for row in b["table_data"]]
    raise RuntimeError("no peptide_table in pepmatch result")


def classify(rows, peptides):
    """rows: pepmatch result rows (only peptides WITH a match are returned).
    Peptides absent from the result have no human match within 1 mismatch.
    similarity_class is a pipeline-defined interpretation, not a validated
    clinical risk taxonomy."""
    out = {p: {"human_match": "none", "min_mismatches": None,
               "similarity_class": "no_match"} for p in peptides}
    for r in rows:
        pep = r["peptide"]
        mm = int(r["mismatches"])
        prot = r.get("protein_name")
        if mm == 0:
            out[pep] = {"human_match": "yes", "min_mismatches": 0,
                        "similarity_class": "exact_match", "best_match_protein": prot}
        else:
            out[pep] = {"human_match": "yes", "min_mismatches": mm,
                        "similarity_class": "one_mismatch", "best_match_protein": prot}
    return out


def main():
    for cls in (1, 2):
        peps = unique_strong_peptides(cls)
        print(f"class {cls}: {len(peps)} unique strong peptides")
        resp = submit(peps)
        pid = resp["pipeline_id"]
        rid = poll(pid)
        rows = fetch_table(rid)
        out = classify(rows, peps)
        out["_meta"] = {
            "tool": "PEPMatch via IEDB next-generation API",
            "version_served": "0.9.5",
            "retrieved": utcnow(),
            "proteome": "Human",
            "mismatch": 1,
            "best_match": True,
            "terminology": "similarity_class is a pipeline-defined sequence-match label (exact_match / one_mismatch / no_match), not a validated clinical risk taxonomy.",
        }
        path = os.path.join(RES, f"human_similarity_class{cls}.json")
        with open(path, "w") as f:
            json.dump(out, f, indent=2)
        print(f"wrote {os.path.basename(path)}")


if __name__ == "__main__":
    main()
