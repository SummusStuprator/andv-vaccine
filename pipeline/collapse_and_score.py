#!/usr/bin/env python3
"""Collapse strong-binder peptides into non-redundant cores, score them, and select the documented analysis set.

The script reads frozen raw prediction CSVs, conservation summaries, the human-proteome sequence-similarity outputs, and the 1000 Genomes HLA typing table. It writes canonical class-I/class-II core tables and selected-set tables.

The scoring formula and selection rule are deterministic and documented in this source file. The sequence-similarity term is a pipeline weighting term, not a clinical safety assessment.
"""
import ast
import glob
import json
import math
import os

import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROT = {"ANDV_Gn": 1, "ANDV_Gc": 2, "ANDV_N": 3}
OFFSET = {"ANDV_Gn": 7, "ANDV_Gc": 477, "ANDV_N": 0}
PANEL = f"{BASE}/data/population/20181129_HLA_types_full_1000_Genomes_Project_panel.txt"

CFG = {
    1: dict(threshold=0.5, n_select=15, el_col="netmhcpan_el_percentile",
            raw_glob="class1_chunk*.csv", safety="human_similarity_class1.json",
            out_cores="cores_class1_scored.csv", out_sel="final_epitopes_class1.csv"),
    2: dict(threshold=2.0, n_select=10, el_col="netmhciipan_el_percentile",
            raw_glob="class2_chunk[0-9].csv", safety="human_similarity_class2.json",
            out_cores="cores_class2_scored.csv", out_sel="final_epitopes_class2.csv"),
}


def load_raw(cls):
    cfg = CFG[cls]
    df = pd.concat([pd.read_csv(f) for f in
                    sorted(glob.glob(f"{BASE}/results/raw_predictions/{cfg['raw_glob']}"))])
    df["protein"] = df["sequence_number"].map({v: k for k, v in PROT.items()})
    if cls == 1:
        df["mhcflurry_processing_score"] = pd.to_numeric(df["mhcflurry_processing_score"],
                                                         errors="coerce")
    else:
        df["score"] = pd.to_numeric(df["score"], errors="coerce")  # CD4Episcore
    return df


def unique_strong(df, cls):
    """One row per unique strong peptide: allele set, best rank, coordinates,
    and the peptide-level sub-score used later."""
    cfg = CFG[cls]
    s = df[df[cfg["el_col"]] < cfg["threshold"]]
    rows = []
    for (protein, peptide), g in s.groupby(["protein", "peptide"]):
        rows.append(dict(
            protein=protein, peptide=peptide,
            start=int(g["start"].iloc[0]), end=int(g["end"].iloc[0]),
            alleles=set(g["allele"]), best_rank=float(g[cfg["el_col"]].min()),
            subscore=(float(g["mhcflurry_processing_score"].iloc[0]) if cls == 1
                      else float(g["score"].iloc[0]))))
    return rows


def overlap_frac(a, b):
    ov = min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1
    shorter = min(a["end"] - a["start"] + 1, b["end"] - b["start"] + 1)
    return ov / shorter


def cluster(peps):
    """Greedy seed-absorption clustering (see module docstring)."""
    peps = sorted(peps, key=lambda p: (-len(p["alleles"]), p["best_rank"], p["start"]))
    assigned, cores = set(), []
    for i, seed in enumerate(peps):
        if i in assigned:
            continue
        members = [seed]
        assigned.add(i)
        for j, cand in enumerate(peps):
            if j in assigned:
                continue
            if overlap_frac(seed, cand) >= 0.60:
                members.append(cand)
                assigned.add(j)
        cores.append(dict(protein=seed["protein"], representative=seed, members=members))
    return cores


def cons_region(cons, protein, start, end):
    off = OFFSET[protein]
    key = "n" if protein == "ANDV_N" else "gpc"
    stats = cons[key]
    lo, hi = start + off, end + off
    rows = [s for s in stats if lo <= s["ref_pos"] <= hi]
    if not rows:
        return None
    def agg(f, e, n):
        return {"mean_frac_identical": round(sum(r[f] for r in rows) / len(rows), 4),
                "min_frac_identical": round(min(r[f] for r in rows), 4),
                "mean_entropy": round(sum(r[e] for r in rows) / len(rows), 4),
                "max_entropy": round(max(r[e] for r in rows), 4),
                "min_covering": min(r[n] for r in rows),
                "n_positions": len(rows)}
    return {"andv": agg("frac_identical_andv", "entropy_andv", "n_covering_andv"),
            "all": agg("frac_identical_all", "entropy_all", "n_covering_all")}


def parse_typed(locus, raw):
    v = str(raw)
    if not v or v in ("None", "nan") or "," in v:
        return None
    a = v.split("/")[0]
    parts = a.split(":")
    if len(parts) > 2:
        a = ":".join(parts[:2])
    return a if a.startswith("HLA") else f"HLA-{a}"


def load_samples():
    pan = pd.read_csv(PANEL, sep="\t")
    samples = []
    for _, r in pan.iterrows():
        c1, c2 = set(), set()
        for loc in ("A", "B", "C"):
            for h in (1, 2):
                a = parse_typed(loc, r[f"HLA-{loc} {h}"])
                if a:
                    c1.add(a)
        for h in (1, 2):
            a = parse_typed("DRB1", r[f"HLA-DRB1 {h}"])
            if a:
                c2.add(a)
            a = parse_typed("DQB1", r[f"HLA-DQB1 {h}"])
            if a:
                c2.add(a)
        samples.append((c1, c2))
    return samples


def c2_credit(alleles):
    """Typed-loci credit: DRB1 direct + DQB1 beta chain of DQ heterodimers."""
    out = set()
    for a in alleles:
        if a.startswith("HLA-DRB1*"):
            out.add(a)
        elif a.startswith("HLA-DQA1*") and "/DQB1*" in a:
            out.add("HLA-DQB1*" + a.split("/DQB1*")[1])
    return out


def build_class(cls):
    cfg = CFG[cls]
    df = load_raw(cls)
    peps = unique_strong(df, cls)
    cons = json.load(open(f"{BASE}/results/conservation_per_position.json"))
    safety = json.load(open(f"{BASE}/results/{cfg['safety']}"))

    cores = []
    for protein in ("ANDV_Gn", "ANDV_Gc", "ANDV_N"):
        cores.extend(cluster([p for p in peps if p["protein"] == protein]))

    # core-level features
    recs = []
    for c in cores:
        mem = c["members"]
        alleles = set().union(*[m["alleles"] for m in mem])
        best_rank = min(m["best_rank"] for m in mem)
        rep = c["representative"]
        region = (min(m["start"] for m in mem), max(m["end"] for m in mem))
        cs = cons_region(cons, c["protein"], *region)
        sims = [safety.get(m["peptide"], {}).get("similarity_class", "no_match") for m in mem]
        recs.append(dict(
            protein=c["protein"], representative_peptide=rep["peptide"],
            rep_start=rep["start"], rep_end=rep["end"],
            region_start=region[0], region_end=region[1],
            n_member_peptides=len(mem), n_alleles_strong=len(alleles),
            alleles_strong=sorted(alleles), best_rank=best_rank,
            cons_ANDV=cs["andv"], cons_ANDV_mean_ident=cs["andv"]["mean_frac_identical"],
            cons_ANDV_min_ident=cs["andv"]["min_frac_identical"],
            cons_all=cs["all"], cons_all_mean_ident=cs["all"]["mean_frac_identical"],
            cons_all_min_ident=cs["all"]["min_frac_identical"],
            member_subscore_max=max(m["subscore"] for m in mem if not math.isnan(m["subscore"]))
                               if any(not math.isnan(m["subscore"]) for m in mem) else float("nan"),
            rep_subscore=rep["subscore"],
            has_exact_human_match=any(s == "exact_match" for s in sims),
            has_onemismatch_human_match=any(s == "one_mismatch" for s in sims),
        ))
    cores_df = pd.DataFrame(recs)

    # scores
    thr = cfg["threshold"]
    max_alleles = cores_df["n_alleles_strong"].max()
    cores_df["binding_score"] = (1 - cores_df["best_rank"] / thr).clip(lower=0)
    cores_df["breadth_score"] = cores_df["n_alleles_strong"] / max_alleles
    cores_df["conservation_score"] = cores_df["cons_ANDV_mean_ident"]
    cores_df["nw_conservation_score"] = cores_df["cons_all_mean_ident"]
    if cls == 1:
        v = cores_df["member_subscore_max"]
        cores_df["processing_score"] = (v - v.min()) / (v.max() - v.min())
    else:
        v = cores_df["rep_subscore"]
        cores_df["immunogenicity_score"] = (v - v.min()) / (v.max() - v.min())
    sub = "processing_score" if cls == 1 else "immunogenicity_score"
    cores_df["similarity_penalty"] = cores_df["has_onemismatch_human_match"].map({True: 0.4, False: 0.0})
    cores_df["total_score"] = (0.30 * cores_df["binding_score"]
                               + 0.25 * cores_df["breadth_score"]
                               + 0.20 * cores_df["conservation_score"]
                               + 0.10 * cores_df[sub]
                               + 0.10 * cores_df["nw_conservation_score"]
                               + 0.05 * (1 - cores_df["similarity_penalty"]))
    cores_df = cores_df.sort_values("total_score", ascending=False).reset_index(drop=True)
    cores_df.insert(0, "core_id", [f"C{cls}-{r['protein']}-{r['rep_start']:03d}"
                                   for _, r in cores_df.iterrows()])
    return cores_df


def select(cores_df, cls):
    """Select the fixed analysis set by composite score after exclusions."""
    cfg = CFG[cls]
    pool = cores_df[~cores_df["has_exact_human_match"]].copy()
    sig = pool["protein"].isin(["ANDV_Gn", "ANDV_Gc"]) & (pool["region_start"] <= 23)
    pool = pool[~sig]
    pool = pool[pool["total_score"].notna()]
    return pool.sort_values("total_score", ascending=False).head(cfg["n_select"])["core_id"].tolist()


def main():
    for cls in (1, 2):
        cfg = CFG[cls]
        cores_df = build_class(cls)
        picked = select(cores_df, cls)
        sel_df = cores_df[cores_df["core_id"].isin(picked)].copy()
        cores_df.to_csv(os.path.join(BASE, "results", cfg["out_cores"]), index=False)
        sel_df.to_csv(os.path.join(BASE, "results", cfg["out_sel"]), index=False)
        print(f"class {cls}: {len(cores_df)} cores; {len(sel_df)} selected")


if __name__ == "__main__":
    main()
