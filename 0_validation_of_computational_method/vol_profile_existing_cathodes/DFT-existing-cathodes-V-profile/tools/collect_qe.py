"""Rebuild QE/ in the clean copy from QE_*.bak + QE/ of the original (read-only on original)."""
import sys, shutil
from pathlib import Path

src_sys = Path(sys.argv[1]).resolve()   # original DFT-<sys>
dst_sys = Path(sys.argv[2]).resolve()   # clean DFT-<sys>

# chronological order: backups by timestamp, then current QE/ last
sources = sorted(src_sys.glob("QE_*.bak")) + [src_sys / "QE"]
cands = {}  # config -> list of (source, status), oldest first
for s in sources:
    for c in sorted(p for p in s.iterdir() if p.is_dir()):
        pwos = list(c.glob("*espresso.pwo"))
        if not pwos:
            continue  # empty / never started
        done = any("JOB DONE" in p.read_text(errors="ignore") for p in pwos)
        cands.setdefault(c.name, []).append((c, done))

qe = dst_sys / "QE"
qe.mkdir(exist_ok=True)
lines = ["config\tstatus\tsource (original folder)"]
for name in sorted(cands):
    lst = cands[name]
    done = [c for c, d in lst if d]
    chosen, status = (done[-1], "complete") if done else (lst[-1][0], "INCOMPLETE")
    t = qe / name
    if t.exists():
        shutil.rmtree(t)
    shutil.copytree(chosen, t, ignore=shutil.ignore_patterns("temp.ipynb", ".*.swp"))
    lines.append(f"{name}\t{status}\t{chosen.relative_to(src_sys)}")
(qe / "SOURCES.tsv").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
