#!/usr/bin/env python3
"""Fetch New World hantavirus M-segment (GPC) and S-segment (N) CDS translations from GenBank.

Reproducible acquisition of the conservation panel.

Design (two-pass, order-independent):
  1. For each query organism name, esearch nuccore:
       "<name>"[Organism] AND <segment terms>
     Segment terms (M): "segment M" OR "M segment" OR glycoprotein OR GPC (all fields).
     Segment terms (S): "segment S" OR "S segment" OR nucleocapsid OR nucleoprotein.
     Note: NCBI Taxonomy maps several historical New World hantavirus names (e.g.
     Lechiguanas, Oran) onto overlapping taxonomic nodes, so queries overlap; each
     record is therefore labeled by its OWN esummary organism annotation, not the query.
  2. All candidate CDS translations are pooled, sorted globally by accession,
     deduplicated by SHA-256 of the amino-acid string, length-filtered
     (GPC 450-1250 aa, retaining partial glycoprotein precursors; N 300-480 aa),
     and capped at 25 per annotated species.

Outputs: data/conservation/{gpc,n}_panel.fasta + fetch_manifest.json
(query strings, fetch timestamp, per-accession organism/length/checksum).
"""
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = Path(__file__).resolve().parent.parent / "data" / "conservation"
OUT.mkdir(parents=True, exist_ok=True)

QUERY_NAMES = [
    "Andes virus", "Sin Nombre virus", "Choclo orthohantavirus", "Rio Negro virus",
    "Laguna Negra virus", "Lechiguanas virus", "Juquitiba virus",
    "Maciel virus", "Oran virus", "Bermejo virus",
]
SEGMENTS = {
    "gpc": {"terms": '("segment M"[All Fields] OR "M segment"[All Fields] OR glycoprotein[All Fields] OR GPC[All Fields])',
            "min_len": 450, "max_len": 1250},
    "n": {"terms": '("segment S"[All Fields] OR "S segment"[All Fields] OR nucleocapsid[All Fields] OR nucleoprotein[All Fields])',
          "min_len": 300, "max_len": 480},
}
CAP_PER_SPECIES = 25


def get(url, **kw):
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=kw.pop("timeout", 60), **kw)
            r.raise_for_status()
            return r
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 * (attempt + 1))


def post(url, **kw):
    for attempt in range(4):
        try:
            r = requests.post(url, timeout=kw.pop("timeout", 300), **kw)
            r.raise_for_status()
            return r
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 * (attempt + 1))


def esearch(query, retmax=1000):
    r = get(f"{EUTILS}/esearch.fcgi", params={
        "db": "nuccore", "term": query, "retmax": retmax, "retmode": "json"})
    return r.json()["esearchresult"]["idlist"]


def efetch_cds_aa(ids):
    """Returns list of (accession, uid, seq)."""
    recs = []
    for i in range(0, len(ids), 200):
        r = post(f"{EUTILS}/efetch.fcgi", data={
            "db": "nuccore", "id": ",".join(ids[i:i + 200]),
            "rettype": "fasta_cds_aa", "retmode": "text"})
        text = r.text
        header, seq = None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if header:
                    recs.append((header, "".join(seq)))
                header, seq = line[1:], []
            else:
                seq.append(line.strip())
        if header:
            recs.append((header, "".join(seq)))
        time.sleep(0.4)
    return recs


def main():
    manifest = {"fetched_at": datetime.now(timezone.utc).isoformat(),
                "source": "NCBI nuccore via E-utilities",
                "query_names": QUERY_NAMES,
                "taxonomy_note": ("NCBI Taxonomy maps several historical New World hantavirus "
                                  "names onto overlapping nodes; records are labeled by their own "
                                  "GenBank organism annotation, not by the query that returned them."),
                "cap_per_species": CAP_PER_SPECIES,
                "length_filters": {k: [v["min_len"], v["max_len"]] for k, v in SEGMENTS.items()},
                "segments": {}, "queries": {}}
    for seg, cfg in SEGMENTS.items():
        # pass 1: collect all candidate ids across query names
        all_ids = {}
        for name in QUERY_NAMES:
            q = f'"{name}"[Organism] AND {cfg["terms"]}'
            manifest["queries"][f"{seg}|{name}"] = q
            try:
                ids = esearch(q)
            except Exception as e:
                print(f"esearch failed {name} {seg}: {e}")
                continue
            for uid in ids:
                all_ids[uid] = name
            time.sleep(0.4)
        ids = sorted(all_ids, key=int)
        recs = efetch_cds_aa(ids)
        # pass 2: global sort by accession, length filter, global dedupe
        seen_hashes, per_species, panel = set(), {}, []
        for header, seq in sorted(recs, key=lambda x: x[0]):
            if not (cfg["min_len"] <= len(seq) <= cfg["max_len"]):
                continue
            h = hashlib.sha256(seq.encode()).hexdigest()
            if h in seen_hashes:
                continue
            acc = header.split()[0]
            acc_base = acc.split("|")[1].split("_prot")[0] if "|" in acc else acc
            seen_hashes.add(h)
            panel.append([acc_base, acc, seq, h])
        # resolve each record's own organism annotation via accession -> esummary
        acc_list = [p[0] for p in panel]
        acc_org = {}
        for i in range(0, len(acc_list), 50):
            r = get(f"{EUTILS}/esearch.fcgi", params={
                "db": "nuccore", "term": " OR ".join(f"{a}[Accession]" for a in acc_list[i:i + 50]),
                "retmode": "json", "retmax": 1000})
            uids = r.json()["esearchresult"]["idlist"]
            time.sleep(0.4)
            if uids:
                r2 = get(f"{EUTILS}/esummary.fcgi", params={
                    "db": "nuccore", "id": ",".join(uids), "retmode": "json"})
                d = r2.json()["result"]
                for uid in uids:
                    if uid in d:
                        acc_org[d[uid]["accessionversion"]] = d[uid].get("organism", "unknown")
                time.sleep(0.4)
        final = []
        for acc_base, acc, seq, h in panel:
            organism = acc_org.get(acc_base, "unknown")
            n_sp = per_species.get(organism, 0)
            if n_sp >= CAP_PER_SPECIES:
                continue
            per_species[organism] = n_sp + 1
            final.append((acc_base, organism, seq, h))
        manifest["segments"][seg] = {
            "n_sequences": len(final),
            "per_species": per_species,
            "accessions": [{"accession": a, "organism": o, "length_aa": len(sq), "sha256": h}
                           for a, o, sq, h in final]}
        with open(OUT / f"{seg}_panel.fasta", "w") as f:
            for a, o, sq, h in final:
                f.write(f">{a}|{o.replace(' ', '_')}\n{sq}\n")
        print(f"{seg}: {len(final)} sequences; per annotated species: {per_species}")
    json.dump(manifest, open(OUT / "fetch_manifest.json", "w"), indent=2)


if __name__ == "__main__":
    main()
