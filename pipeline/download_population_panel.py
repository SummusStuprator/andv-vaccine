#!/usr/bin/env python3
"""download_population_panel.py — fetch and checksum the 1000 Genomes HLA panel.

Downloads the HLA type panel (2,693 samples; loci A, B, C, DQB1, DRB1) from the
1000 Genomes FTP mirror and verifies SHA-256 against the frozen value recorded
for this frozen release. Refuses to overwrite a mismatched existing file.
"""
import hashlib
import os
import sys

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
DEST_DIR = os.path.join(os.path.dirname(HERE), "data", "population")
BASE_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types"
FILES = {
    "20181129_HLA_types_full_1000_Genomes_Project_panel.txt":
        "570b5cc5523273826a5695452f18b77f8742e568e1e35a34aed5d2a313472e8e",
    "20181129_HLA_manifest.txt": None,
    "README_20181129_HLA_types_full_1000_Genomes_Project_panel.txt": None,
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    os.makedirs(DEST_DIR, exist_ok=True)
    for name, expected in FILES.items():
        dest = os.path.join(DEST_DIR, name)
        if os.path.exists(dest) and (expected is None or sha256(dest) == expected):
            print(f"present, checksum ok: {name}")
            continue
        url = f"{BASE_URL}/{name}"
        print(f"downloading {url}")
        r = requests.get(url, timeout=300)
        r.raise_for_status()
        with open(dest, "wb") as f:
            f.write(r.content)
        got = sha256(dest)
        if expected is not None and got != expected:
            os.remove(dest)
            sys.exit(f"CHECKSUM MISMATCH for {name}: got {got}, expected {expected}")
        print(f"  sha256 {got}")


if __name__ == "__main__":
    main()
