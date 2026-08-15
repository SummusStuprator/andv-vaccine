#!/usr/bin/env python3
"""Post-process archived BepiPred and DiscoTope per-residue outputs.

This script starts from frozen external predictor tables in `results/` and produces a compact B-cell prediction summary. It does not establish neutralization, and it does not recreate the unavailable upstream structure-to-DiscoTope stage.
"""
import json
import os
import re

import pandas as pd
from Bio import SeqIO

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")

MIN_LINEAR_LEN = 5
RSA_MIN = 0.2
PLDDT_MIN = 70
PATCH_GAP = 2
OFFSET = {"ANDV_Gn": 7, "ANDV_Gc": 477, "ANDV_N": 0}


def linear_regions(df):
    e = df["bepipred_assignment"].values
    regions, start = [], None
    for i, a in enumerate(e):
        if a == "E" and start is None:
            start = i
        elif a != "E" and start is not None:
            regions.append((start + 1, i))
            start = None
    if start is not None:
        regions.append((start + 1, len(e)))
    return [{"start": s, "end": e, "length": e - s + 1}
            for s, e in regions if e - s + 1 >= MIN_LINEAR_LEN]


def patches(df):
    hc = df[(df["assignment"] == "E") & (df["rsa"] > RSA_MIN)
            & (df["discotope.plddts"] > PLDDT_MIN)]
    nums = sorted(hc["res_num"])
    out, start, prev = [], None, None
    for n in nums:
        if start is None:
            start = prev = n
        elif n - prev <= PATCH_GAP + 1:
            prev = n
        else:
            out.append({"start": start, "end": prev, "n_residues": prev - start + 1})
            start = prev = n
    if start is not None:
        out.append({"start": start, "end": prev, "n_residues": prev - start + 1})
    return len(hc), out


def sequons(seq):
    return [m.start() + 1 for m in re.finditer(r"N[^P][ST]", seq)]


def region_cons(cons, key, lo, hi):
    rows = [s for s in cons[key] if lo <= s["ref_pos"] <= hi]
    return round(sum(r["frac_identical_andv"] for r in rows) / len(rows), 4)


def main():
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(
        os.path.join(ROOT, "data", "canonical_sequences.fasta"), "fasta")}
    cons = json.load(open(os.path.join(RES, "conservation_per_position.json")))

    out = {"terminology_note": ("BepiPred-3.0/DiscoTope-3.0 predict B-cell epitopes, "
                                "not neutralization."),
           "bepipred_linear_regions": {}, "discotope": {}, "nglyco_sequons": {},
           "modules": {}}
    for p in ("ANDV_Gn", "ANDV_Gc", "ANDV_N"):
        bp = pd.read_csv(os.path.join(RES, f"bepipred_{p}.csv"))
        out["bepipred_linear_regions"][p] = linear_regions(bp)
        out["nglyco_sequons"][p] = sequons(seqs[p])
    for p in ("ANDV_Gn", "ANDV_Gc"):
        dt = pd.read_csv(os.path.join(RES, f"discotope_{p}.csv"))
        n_hc, pchs = patches(dt)
        out["discotope"][p] = {
            "n_high_confidence_residues": n_hc,
            "plddt_column_note": ("discotope.plddts is uniformly 100.0 in the shipped "
                                  "CSV; the pLDDT>70 filter is vacuous on these data."),
            "sequence_contiguity_patches": pchs}

    # module profiles
    gn_lo, gn_hi = 272, 287
    gc_seq = seqs["ANDV_Gc"]
    gc_lo = gc_seq.find("VTGFNQIDSDKVYDD") + 1
    gc_hi = gc_lo + 14
    dtg = pd.read_csv(os.path.join(RES, "discotope_ANDV_Gn.csv"))
    sub = dtg[(dtg["res_num"] >= gn_lo) & (dtg["res_num"] <= gn_hi)]
    out["modules"]["Gn_272_287"] = {
        "sequence": seqs["ANDV_Gn"][gn_lo - 1:gn_hi],
        "bepipred_region": True,
        "discotope_E_fraction": round((sub["assignment"] == "E").mean(), 3),
        "surface_exposed_fraction_rsa>0.2": round((sub["rsa"] > RSA_MIN).mean(), 3),
        "cons_ANDV_mean_ident": region_cons(cons, "gpc", gn_lo + 7, gn_hi + 7),
        "glycan_shielded": any(gn_lo <= s <= gn_hi for s in out["nglyco_sequons"]["ANDV_Gn"])}
    dtc = pd.read_csv(os.path.join(RES, "discotope_ANDV_Gc.csv"))
    subc = dtc[(dtc["res_num"] >= gc_lo) & (dtc["res_num"] <= gc_hi)]
    out["modules"]["Gc_591_605"] = {
        "sequence": "VTGFNQIDSDKVYDD",
        "stored_start": gc_lo, "stored_end": gc_hi,
        "bepipred_region": False,
        "discotope_E_fraction": round((subc["assignment"] == "E").mean(), 3),
        "surface_exposed_fraction_rsa>0.2": round((subc["rsa"] > RSA_MIN).mean(), 3),
        "cons_ANDV_mean_ident": region_cons(cons, "gpc", gc_lo + 477, gc_hi + 477),
        "glycan_shielded": any(gc_lo <= s <= gc_hi for s in out["nglyco_sequons"]["ANDV_Gc"]),
        "note": "same sequence as class-II core C2-ANDV_Gc-591 representative peptide"}

    with open(os.path.join(RES, "bcell_module.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: (len(v) if isinstance(v, list) else v)
                      for k, v in out["bepipred_linear_regions"].items()}, indent=1))
    print("Gn discotope HC:", out["discotope"]["ANDV_Gn"]["n_high_confidence_residues"],
          "| Gc:", out["discotope"]["ANDV_Gc"]["n_high_confidence_residues"])
    print("sequons:", out["nglyco_sequons"])
    print("modules:", json.dumps(out["modules"], indent=1))


if __name__ == "__main__":
    main()
