"""Build a clean copy of DFT-existing-cathodes-V-profile. Original is only read."""
import os, sys, shutil, tarfile, json
from pathlib import Path

SRC = Path(sys.argv[1]).resolve()
DST = Path(sys.argv[2]).resolve()
assert not DST.exists(), f"{DST} already exists"

KEEP_TOP_FILES = {"settings.yaml", "full_simulation_settings.yaml", "job.sh", "simulation.log"}
KEEP_TOP_DIRS = {"QE", "plots"}
KEEP_VP = {"voltage.csv", "hull.csv", "thermo.csv", "trj.xyz", "states"}
ANALYSIS = {"vis.ipynb", "temp.ipynb", "pbe_settings.yaml",
            "VP/plt_hull.py", "VP/formation_energy_plot.png", "QE/temp.ipynb"}


def is_junk(name):
    return name == ".DS_Store" or name.startswith("._") or name.endswith(".swp")


def classify(rel):
    """rel: path relative to system dir (posix str) -> keep | analysis | archive | junk"""
    parts = rel.split("/")
    if any(is_junk(p) for p in parts):
        return "junk"
    if rel in ANALYSIS:
        return "analysis"
    top = parts[0]
    if "bak" in top:
        return "archive"
    if len(parts) == 1:
        if top in KEEP_TOP_FILES or top.endswith(".cif"):
            return "keep"
        return "archive"
    if top in KEEP_TOP_DIRS:
        return "keep"
    if top == "VP":
        return "keep" if parts[1] in KEEP_VP else "archive"
    return "archive"


DST.mkdir()
(DST / "archive").mkdir()
report = {}
for sysdir in sorted(SRC.glob("DFT-*")):
    out = DST / sysdir.name
    groups = {"keep": [], "analysis": [], "archive": [], "junk": []}
    for p in sorted(sysdir.rglob("*")):
        if p.is_file() or p.is_symlink():
            rel = p.relative_to(sysdir).as_posix()
            groups[classify(rel)].append(rel)
    for rel in groups["keep"]:
        t = out / rel
        t.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(sysdir / rel, t)
    for rel in groups["analysis"]:
        t = out / "analysis" / Path(rel).name
        t.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(sysdir / rel, t)
    if groups["archive"]:
        with tarfile.open(DST / "archive" / f"{sysdir.name}_old_runs.tar.gz", "w:gz") as tf:
            for rel in groups["archive"]:
                tf.add(sysdir / rel, arcname=f"{sysdir.name}/{rel}", recursive=False)
    report[sysdir.name] = {k: len(v) for k, v in groups.items()}
    report[sysdir.name]["junk_files"] = groups["junk"]
print(json.dumps(report, indent=1))
