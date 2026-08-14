#!/usr/bin/env python3
"""Conservation analysis: align the GenBank panel with MAFFT and score per-position
and per-epitope-core conservation relative to the native ANDV reference.

Inputs:
  data/conservation/gpc_panel.fasta, data/conservation/n_panel.fasta   (from fetch_conservation_sequences.py)
  data/native_genbank_reference.fasta                                  (ANDV_GPC_native 1138 aa, ANDV_N_native 428 aa)
  results/cores_class1_scored.csv, results/cores_class2_scored.csv     (core coordinates on STORED sequences)

Coordinate mapping (stored -> native):
  stored Gn = 23-aa engineered signal + GPC[31-500]   => native GPC pos = stored pos + 7
  stored Gc = 23-aa engineered signal + GPC[501-1138] => native GPC pos = stored pos + 477
  stored N  = native N (identical)

Per-position statistics (gaps excluded; positions with < MIN_COVER covering sequences
are reported with n_covering and flagged low_coverage=True):
  frac_identical to the ANDV reference residue, Shannon entropy (base 2),
  computed twice: over ANDV-annotated sequences only, and over the full panel.

Outputs:
  results/conservation_per_position.json
  results/cores_class{1,2}_conservation.csv  (cores + recomputed conservation columns)
"""
import csv
import json
import math
import os
import subprocess
import tempfile
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
RES = os.path.join(ROOT, "results")
MIN_COVER = 3

# stored-sequence coordinate offset -> native reference position (1-based, inclusive)
OFFSET = {"ANDV_Gn": 7, "ANDV_Gc": 477, "ANDV_N": 0}
REFKEY = {"ANDV_Gn": "gpc", "ANDV_Gc": "gpc", "ANDV_N": "n"}


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif line:
            seqs[name] += line
    return seqs


def mafft_align(seqs):
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as f:
        for k, v in seqs.items():
            f.write(f">{k}\n{v}\n")
        inf = f.name
    out = subprocess.run(["mafft", "--auto", inf], capture_output=True, text=True, check=True)
    os.unlink(inf)
    return read_fasta_lines(out.stdout.splitlines())


def read_fasta_lines(lines):
    seqs, name = {}, None
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = ""
        elif line:
            seqs[name] += line
    return seqs


def per_position_stats(aln, ref_name):
    """Return list over reference ungapped positions: dict with frac_identical/entropy
    for ANDV-annotated subset and full panel."""
    ref = aln[ref_name]
    names = [n for n in aln if n != ref_name]
    andv_names = [n for n in names if "Andes" in n or "andesense" in n]
    stats = []
    for i, aa in enumerate(ref):
        if aa == "-":
            continue
        def col_stats(subset):
            chars = [aln[n][i] for n in subset if aln[n][i] != "-"]
            n_cov = len(chars)
            if n_cov == 0:
                return 0.0, 0.0, 0
            frac = sum(1 for c in chars if c == aa) / n_cov
            counts = Counter(chars)
            ent = -sum((c / n_cov) * math.log2(c / n_cov) for c in counts.values())
            return frac, ent, n_cov
        f_all, e_all, n_all = col_stats(names)
        f_andv, e_andv, n_andv = col_stats(andv_names)
        stats.append({
            "ref_pos": len(stats) + 1, "ref_aa": aa,
            "frac_identical_all": round(f_all, 4), "entropy_all": round(e_all, 4),
            "n_covering_all": n_all, "low_coverage_all": n_all < MIN_COVER,
            "frac_identical_andv": round(f_andv, 4), "entropy_andv": round(e_andv, 4),
            "n_covering_andv": n_andv, "low_coverage_andv": n_andv < MIN_COVER,
        })
    return stats


def core_stats(stats, start_ref, end_ref):
    """Aggregate per-position stats over a core's reference coordinate span."""
    rows = [s for s in stats if start_ref <= s["ref_pos"] <= end_ref]
    if not rows:
        return None
    def agg(fkey, ekey, nkey):
        covered = [r for r in rows if not r[fkey.replace("frac_identical", "low_coverage").replace("identical", "coverage")] ] if False else rows
        fr = [r[fkey] for r in rows]
        en = [r[ekey] for r in rows]
        nc = [r[nkey] for r in rows]
        return {
            "mean_frac_identical": round(sum(fr) / len(fr), 4),
            "min_frac_identical": round(min(fr), 4),
            "mean_entropy": round(sum(en) / len(en), 4),
            "max_entropy": round(max(en), 4),
            "min_covering": min(nc),
            "n_positions": len(rows),
        }
    return {"andv": agg("frac_identical_andv", "entropy_andv", "n_covering_andv"),
            "all": agg("frac_identical_all", "entropy_all", "n_covering_all")}


def main():
    native = read_fasta(os.path.join(DATA, "native_genbank_reference.fasta"))
    panels = {"gpc": read_fasta(os.path.join(DATA, "conservation", "gpc_panel.fasta")),
              "n": read_fasta(os.path.join(DATA, "conservation", "n_panel.fasta"))}
    refs = {"gpc": ("ANDV_reference_GPC|MN258217.1", native["ANDV_GPC_native"]),
            "n": ("ANDV_reference_N|MN258242.1", native["ANDV_N_native"])}
    all_stats = {}
    for seg in ["gpc", "n"]:
        ref_name, ref_seq = refs[seg]
        seqs = dict(panels[seg])
        seqs[ref_name] = ref_seq
        print(f"aligning {seg}: {len(seqs)} sequences ...")
        aln = mafft_align(seqs)
        all_stats[seg] = per_position_stats(aln, ref_name)
        print(f"  {len(all_stats[seg])} reference positions scored")
    json.dump(all_stats, open(os.path.join(RES, "conservation_per_position.json"), "w"))

    for cls in [1, 2]:
        inp = os.path.join(RES, f"cores_class{cls}_scored.csv")
        outp = os.path.join(RES, f"cores_class{cls}_conservation.csv")
        with open(inp) as f, open(outp, "w", newline="") as g:
            rows = list(csv.DictReader(f))
            w = csv.writer(g)
            w.writerow(["core_id", "protein", "representative_peptide", "rep_start", "rep_end",
                        "ref_start", "ref_end",
                        "cons_andv_mean", "cons_andv_min", "cons_andv_mincover",
                        "cons_all_mean", "cons_all_min", "cons_all_mincover"])
            for r in rows:
                off = OFFSET[r["protein"]]
                seg = REFKEY[r["protein"]]
                rs, re_ = int(r["rep_start"]) + off, int(r["rep_end"]) + off
                cs = core_stats(all_stats[seg], rs, re_)
                if cs is None:
                    continue
                w.writerow([r["core_id"], r["protein"], r["representative_peptide"],
                            r["rep_start"], r["rep_end"], rs, re_,
                            cs["andv"]["mean_frac_identical"], cs["andv"]["min_frac_identical"], cs["andv"]["min_covering"],
                            cs["all"]["mean_frac_identical"], cs["all"]["min_frac_identical"], cs["all"]["min_covering"]])
        print(f"wrote {outp} ({len(rows)} cores)")


if __name__ == "__main__":
    main()
