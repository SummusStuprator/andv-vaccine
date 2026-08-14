#!/usr/bin/env python3
"""Build MANIFEST.sha256 from stable Git-tracked release files."""
from __future__ import annotations
import hashlib
import subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'MANIFEST.sha256'

def tracked_files():
    try:
        raw=subprocess.check_output(['git','-C',str(ROOT),'ls-files','-z'])
        paths=[ROOT/p.decode() for p in raw.split(b'\0') if p]
        if paths:
            return paths
    except Exception:
        pass
    return [p for p in ROOT.rglob('*') if p.is_file() and '.git' not in p.parts and '__pycache__' not in p.parts and p.suffix != '.pyc']

rows=[]
for p in sorted(tracked_files()):
    if p == OUT or not p.is_file():
        continue
    rel=p.relative_to(ROOT).as_posix()
    rows.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {rel}")
OUT.write_text('\n'.join(rows)+'\n',encoding='utf-8')
print(f"wrote MANIFEST.sha256: {len(rows)} files")
