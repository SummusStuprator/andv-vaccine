#!/usr/bin/env python3
"""Conservation analysis: align the GenBank panels with MAFFT and score per-position
conservation relative to the native ANDV reference.

Inputs:
  data/conservation/gpc_panel.fasta, data/conservation/n_panel.fasta
  data/native_genbank_reference.fasta

Per-position statistics (gaps excluded; positions with < MIN_COVER covering sequences
are reported with n_covering and flagged low_coverage=True):
  fraction identical to the ANDV reference residue and Shannon entropy (base 2),
  computed for the ANDV-annotated subset and the full panel.

Output:
  results/conservation_per_position.json

Core-region conservation is calculated downstream by collapse_and_score.py from this
per-position file, using each core's full region coordinates.
"""
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
    """Return statistics for each ungapped reference position."""
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
            "ref_pos": len(stats) + 1,
            "ref_aa": aa,
            "frac_identical_all": round(f_all, 4),
            "entropy_all": round(e_all, 4),
            "n_covering_all": n_all,
            "low_coverage_all": n_all < MIN_COVER,
            "frac_identical_andv": round(f_andv, 4),
            "entropy_andv": round(e_andv, 4),
            "n_covering_andv": n_andv,
            "low_coverage_andv": n_andv < MIN_COVER,
        })
    return stats


def main():
    native = read_fasta(os.path.join(DATA, "native_genbank_reference.fasta"))
    panels = {
        "gpc": read_fasta(os.path.join(DATA, "conservation", "gpc_panel.fasta")),
        "n": read_fasta(os.path.join(DATA, "conservation", "n_panel.fasta")),
    }
    refs = {
        "gpc": ("ANDV_reference_GPC|MN258217.1", native["ANDV_GPC_native"]),
        "n": ("ANDV_reference_N|MN258242.1", native["ANDV_N_native"]),
    }
    all_stats = {}
    for seg in ["gpc", "n"]:
        ref_name, ref_seq = refs[seg]
        seqs = dict(panels[seg])
        seqs[ref_name] = ref_seq
        print(f"aligning {seg}: {len(seqs)} sequences ...")
        aln = mafft_align(seqs)
        all_stats[seg] = per_position_stats(aln, ref_name)
        print(f"  {len(all_stats[seg])} reference positions scored")
    with open(os.path.join(RES, "conservation_per_position.json"), "w", encoding="utf-8") as f:
        json.dump(all_stats, f)
    print("wrote results/conservation_per_position.json")


if __name__ == "__main__":
    main()
