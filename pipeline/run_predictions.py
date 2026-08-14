#!/usr/bin/env python3
"""Regenerate raw IEDB next-gen API epitope predictions for the ANDV vaccine pipeline.

Submits chunked prediction jobs (class I: NetMHCpan EL+BA, MHCflurry, IEDB
immunogenicity; class II: NetMHCIIpan EL+BA, CD4Episcore) on the canonical
stored sequences, polls to completion, and saves per-chunk peptide tables as
CSV under results/raw_predictions/. Resumable: chunks with an existing
non-empty CSV are skipped. All API responses (including any version metadata)
are archived for provenance.

Usage:
    python run_predictions.py --class 1 [--test]   # submit + retrieve class I
    python run_predictions.py --class 2 [--test]   # submit + retrieve class II
"""
import argparse
import csv
import json
import os
import sys
import time
from datetime import datetime, timezone

import requests

BASE = "https://api-nextgen-tools.iedb.org/api/v1"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "raw_predictions")

CLASS1_PREDICTORS = [
    {"type": "binding", "method": "netmhcpan_el"},
    {"type": "binding", "method": "netmhcpan_ba"},
    {"type": "binding", "method": "mhcflurry"},
    {"type": "immunogenicity", "method": "immunogenicity"},
]
CLASS2_PREDICTORS = [
    {"type": "binding", "method": "netmhciipan_el"},
    {"type": "binding", "method": "netmhciipan_ba"},
    {"type": "immunogenicity", "method": "cd4episcore"},
]


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = ""
            elif line:
                seqs[name] += line
    return seqs


def load_alleles(cls):
    panel = json.load(open(os.path.join(ROOT, "data", "hla_allele_panel.json")))["hla_allele_panel"]
    if cls == 1:
        return [a["allele"] for a in panel["class_I"]["alleles"]]
    return list(panel["class_II"]["alleles"])


def submit(title, tool_group, sequence_text, alleles, length_range, predictors):
    payload = {
        "pipeline_title": title,
        "run_stage_range": [1, 1],
        "stages": [{
            "stage_number": 1,
            "tool_group": tool_group,
            "input_sequence_text": sequence_text,
            "input_parameters": {
                "alleles": alleles,
                "peptide_length_range": length_range,
                "predictors": predictors,
            },
        }],
    }
    r = requests.post(f"{BASE}/pipeline", json=payload, timeout=120)
    r.raise_for_status()
    return r.json()


def stage_statuses(pipeline_id):
    r = requests.get(f"{BASE}/pipeline/{pipeline_id}", timeout=60)
    r.raise_for_status()
    stages = r.json().get("stages", [])
    out = []
    for s in stages:
        st = s.get("stage_status")
        if st is None and s.get("stage_url"):
            try:
                st = requests.get(s["stage_url"], timeout=60).json().get("stage_status")
            except Exception:
                st = None
        out.append(st)
    return out


def fetch_peptide_table(result_id):
    r = requests.get(f"{BASE}/results/{result_id}", timeout=300)
    r.raise_for_status()
    d = r.json()
    # archive the raw response envelope (minus bulky table) for provenance
    meta = {"fetched_at": utcnow(), "result_id": result_id,
            "top_level_keys": list(d.keys()),
            "data_keys": list(d.get("data", {}).keys()) if isinstance(d.get("data"), dict) else None,
            "versions": d.get("versions") or d.get("data", {}).get("versions")}
    for b in d["data"]["results"]:
        if b.get("type") == "peptide_table":
            cols = [c["name"] for c in b["table_columns"]]
            return cols, b["table_data"], meta
    return None, None, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", type=int, required=True, choices=[1, 2])
    ap.add_argument("--test", action="store_true", help="single allele, ANDV_N only")
    ap.add_argument("--poll-interval", type=int, default=60)
    args = ap.parse_args()
    cls = args.cls
    os.makedirs(OUT, exist_ok=True)

    seqs = read_fasta(os.path.join(ROOT, "data", "canonical_sequences.fasta"))
    seqs = {k: v for k, v in seqs.items() if k.startswith("ANDV")}
    if args.test:
        seqs = {"ANDV_N": seqs["ANDV_N"]}
    combined = "\n".join(f">{k}\n{v}" for k, v in seqs.items())

    alleles = load_alleles(cls)
    if args.test:
        alleles = alleles[:1]
    chunk = 14 if cls == 1 else 11
    predictors = CLASS1_PREDICTORS if cls == 1 else CLASS2_PREDICTORS
    tool_group = "mhci" if cls == 1 else "mhcii"
    length_range = [8, 11] if cls == 1 else [15, 15]
    tag = f"class{cls}"

    jobs_path = os.path.join(OUT, f"{tag}_jobs.json")
    if os.path.exists(jobs_path):
        jobs = json.load(open(jobs_path))
        print(f"resuming: loaded {len(jobs)} jobs from {jobs_path}")
    else:
        jobs = []
        for i in range(0, len(alleles), chunk):
            ca = alleles[i:i + chunk]
            resp = submit(f"ANDV_publication_{tag}_chunk{i // chunk}", tool_group, combined,
                          ",".join(ca), length_range, predictors)
            rid = resp.get("result_id") or resp.get("id") or resp.get("data", {}).get("result_id")
            jobs.append({"chunk": i // chunk, "alleles": ca, "result_id": rid,
                         "pipeline_id": resp.get("pipeline_id"),
                         "submitted_at": utcnow(), "submit_response": resp})
            print(f"[{utcnow()}] submitted {tag} chunk{i // chunk}: {len(ca)} alleles -> {rid}", flush=True)
            time.sleep(2)
        json.dump(jobs, open(jobs_path, "w"), indent=2)

    pending = list(jobs)
    while pending:
        still = []
        for j in pending:
            outcsv = os.path.join(OUT, f"{tag}_chunk{j['chunk']}.csv")
            if os.path.exists(outcsv) and os.path.getsize(outcsv) > 0:
                continue
            try:
                statuses = stage_statuses(j["pipeline_id"])
            except Exception as e:
                print(f"{tag} chunk{j['chunk']}: status error {e}", flush=True)
                still.append(j)
                continue
            if all(s == "done" for s in statuses):
                try:
                    cols, rows, meta = fetch_peptide_table(j["result_id"])
                except Exception as e:
                    print(f"{tag} chunk{j['chunk']}: fetch error {e}", flush=True)
                    still.append(j)
                    continue
                if cols:
                    with open(outcsv, "w", newline="") as f:
                        w = csv.writer(f)
                        w.writerow(cols)
                        w.writerows(rows)
                    json.dump(meta, open(os.path.join(OUT, f"{tag}_chunk{j['chunk']}_meta.json"), "w"), indent=2)
                    print(f"[{utcnow()}] {tag} chunk{j['chunk']}: saved {len(rows)} rows", flush=True)
                else:
                    print(f"{tag} chunk{j['chunk']}: done but no peptide_table", flush=True)
                    still.append(j)
            else:
                print(f"{tag} chunk{j['chunk']}: status={statuses}", flush=True)
                still.append(j)
        pending = still
        if pending:
            time.sleep(args.poll_interval)
    print("ALL DONE")


if __name__ == "__main__":
    main()
