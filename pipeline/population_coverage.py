#!/usr/bin/env python3
"""Population-coverage analysis from 1000 Genomes per-sample HLA genotypes.

Coverage is evaluated for the representative peptides in the fixed selected set.
The strong-binding allele set for each representative peptide is reconstructed
from the frozen raw prediction tables; cluster-member allele unions are not used
for the coverage estimate.

The Abi-Rached et al. HLA panel types HLA-A, -B, -C, -DRB1 and -DQB1 only.
Class-I coverage therefore uses A/B/C directly. Class-II coverage is reported
under two assumptions: DRB1 only, and DRB1 plus DQB1 beta-chain credit for DQ
heterodimers. Because DQA1 is not typed, the latter is an upper-bound-type
assumption for the DQ contribution. DRB3/4/5 and DP cannot be credited.

A sample is covered for a class when at least one typed allele is among the
strong-binding alleles of at least one selected representative peptide.
Combined coverage requires both class-I and class-II coverage. Bootstrap
confidence intervals resample individuals with replacement.
"""
import csv
import glob
import json
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
RAW = os.path.join(RES, "raw_predictions")
GENO = os.path.join(ROOT, "data", "population",
                    "20181129_HLA_types_full_1000_Genomes_Project_panel.txt")
N_BOOT = 2000
SEED = 20260812


def parse_allele(locus, raw):
    """Normalize a genotype-table entry to HLA-X*NN:NN or None."""
    v = (raw or "").strip()
    if v in ("", "-", "NA", "na", "None", "null"):
        return None
    if "/" in v:
        v = v.split("/")[0]
    if " " in v:
        v = v.split()[0]
    parts = v.split(":")
    if len(parts) < 2 or not parts[0].isdigit():
        return None
    return f"HLA-{locus}*{':'.join(parts[:2])}"


def load_samples():
    samples = []
    with open(GENO, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            c1, drb1, dqb1 = set(), set(), set()
            for loc in ("A", "B", "C"):
                for h in ("1", "2"):
                    a = parse_allele(loc, r[f"HLA-{loc} {h}"])
                    if a:
                        c1.add(a)
            for h in ("1", "2"):
                a = parse_allele("DRB1", r[f"HLA-DRB1 {h}"])
                if a:
                    drb1.add(a)
                a = parse_allele("DQB1", r[f"HLA-DQB1 {h}"])
                if a:
                    dqb1.add(a)
            samples.append({"id": r["Sample ID"], "region": r["Region"],
                            "pop": r["Population"], "c1": c1,
                            "drb1": drb1, "dqb1": dqb1})
    return samples


def load_selected_rep_alleles(cls):
    """Return selected core_id -> strong-binding alleles of its representative peptide."""
    sel_path = os.path.join(RES, f"final_epitopes_class{cls}.csv")
    with open(sel_path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    peptide_to_id = {r["representative_peptide"]: r["core_id"] for r in rows}
    out = {r["core_id"]: set() for r in rows}
    if cls == 1:
        files = sorted(glob.glob(os.path.join(RAW, "class1_chunk*.csv")))
        rank_col, threshold = "netmhcpan_el_percentile", 0.5
    else:
        files = sorted(glob.glob(os.path.join(RAW, "class2_chunk[0-9].csv")))
        rank_col, threshold = "netmhciipan_el_percentile", 2.0
    if not files:
        raise FileNotFoundError("raw prediction CSVs are required; see results/raw_predictions/README.md")
    for path in files:
        with open(path, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                cid = peptide_to_id.get(r["peptide"])
                if cid is not None and float(r[rank_col]) < threshold:
                    out[cid].add(r["allele"])
    missing = [cid for cid, alleles in out.items() if not alleles]
    if missing:
        raise RuntimeError(f"no strong-binding allele rows found for selected representatives: {missing}")
    return out


def c2_typed_locus_map(rep_alleles, mode):
    out = {}
    for cid, alleles in rep_alleles.items():
        drb1, dqb1 = set(), set()
        for a in alleles:
            if a.startswith("HLA-DRB1*"):
                drb1.add(a)
            elif mode == "drb1_dqb1" and a.startswith("HLA-DQA1*") and "/DQB1*" in a:
                dqb1.add("HLA-DQB1*" + a.split("/DQB1*")[1])
        out[cid] = {"drb1": drb1, "dqb1": dqb1}
    return out


def covered_c1(sample, c1_alleles):
    return any(alleles & sample["c1"] for alleles in c1_alleles.values())


def covered_c2(sample, c2_map):
    return any((m["drb1"] & sample["drb1"]) or (m["dqb1"] & sample["dqb1"])
               for m in c2_map.values())


def coverage(samples, c1_alleles, c2_map):
    n = len(samples)
    n1 = sum(covered_c1(s, c1_alleles) for s in samples)
    n2 = sum(covered_c2(s, c2_map) for s in samples)
    nb = sum(covered_c1(s, c1_alleles) and covered_c2(s, c2_map) for s in samples)
    return 100 * n1 / n, 100 * n2 / n, 100 * nb / n


def bootstrap(samples, c1_alleles, c2_map, n_boot=N_BOOT, seed=SEED):
    rng = random.Random(seed)
    n = len(samples)
    flags = [covered_c1(s, c1_alleles) and covered_c2(s, c2_map) for s in samples]
    boots = []
    for _ in range(n_boot):
        hits = sum(flags[rng.randrange(n)] for _ in range(n))
        boots.append(100 * hits / n)
    boots.sort()
    return boots[int(0.025 * n_boot)], boots[int(0.975 * n_boot)]


def main():
    samples = load_samples()
    c1 = load_selected_rep_alleles(1)
    c2 = load_selected_rep_alleles(2)
    result = {
        "model": "per-sample genotype coverage of selected representative peptides",
        "n_samples": len(samples),
        "bootstrap_replicates": N_BOOT,
        "seed": SEED,
        "allele_credit": "representative-peptide strong-binding allele sets from frozen raw predictions",
        "scope_note": ("Class II is restricted to typed loci: DRB1 directly and, in the broader "
                       "variant, the DQB1 beta chain of DQ heterodimers. DQA1 is untyped; "
                       "DRB3/4/5 and DP are not credited."),
        "selected_core_ids": {"class_I": list(c1), "class_II": list(c2)},
        "representative_allele_breadth": {
            "class_I": {cid: len(a) for cid, a in c1.items()},
            "class_II": {cid: len(a) for cid, a in c2.items()},
        },
    }
    for mode in ("drb1_only", "drb1_dqb1"):
        c2_map = c2_typed_locus_map(c2, mode)
        p1, p2, pb = coverage(samples, c1, c2_map)
        lo, hi = bootstrap(samples, c1, c2_map)
        by_region = {}
        for reg in sorted({s["region"] for s in samples}):
            sub = [s for s in samples if s["region"] == reg]
            r1, r2, rb = coverage(sub, c1, c2_map)
            by_region[reg] = {"class_I": round(r1, 2), "class_II": round(r2, 2),
                              "combined": round(rb, 2), "n": len(sub)}
        result[mode] = {
            "class_I_pct": round(p1, 2), "class_II_pct": round(p2, 2),
            "combined_pct": round(pb, 2), "combined_95CI": [round(lo, 2), round(hi, 2)],
            "by_superpopulation_combined": by_region,
        }
        print(f"{mode}: C1={p1:.2f} C2={p2:.2f} combined={pb:.2f} (95% CI {lo:.2f}-{hi:.2f})")
    with open(os.path.join(RES, "population_coverage.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print("wrote population_coverage.json")


if __name__ == "__main__":
    main()
