"""Recompute convex hull and 0 K voltage steps for each DFT-<sys>/ in the clean folder.

Reimplements the pymc zerok-voltage hull/voltage step (formation energy per host, lower convex
hull, V = -[E(n2) - E(n1) - (n2-n1) mu] / ((n2-n1) z)). Checked: fed with VP/states energies it
reproduces VP/voltage.csv for KVPO4F, LiCoO2, MgV2O4, NaCoO2 (for Li2MnO3 and LiFePO4 some
VP/states entries are stale; VP/hull.csv shows the driver used the relaxed values).

  original  VP/voltage.csv as written by the driver
  A         state energy = lowest finished QE relaxation of that composition, all with ONE
            Ry->eV constant (ASE); VP/states entries without any QE output are kept.
            mu as in settings.yaml.
  B         A + mu from a finished bulk-metal QE run with the same constant (Li systems without
            their own finished Li run -> Li2MnO3's Li run)

Usage: python tools/recompute_voltage.py <clean_root>
Writes <clean_root>/recomputed/{<sys>_voltage.csv, <sys>_hull.csv, voltage_comparison.csv}.
"""
import csv, re, sys
from pathlib import Path

RY = 13.605693012183622  # ase.units.Ry (CODATA 2014)
ION_OF = {"KVPO4F": "K", "Li2MnO3": "Li", "LiCoO2": "Li", "LiFePO4": "Li", "MgV2O4": "Mg", "NaCoO2": "Na"}


def n_ion(formula, ion):
    m = re.search(rf"{ion}(\d*)(?=[A-Z]|$)", formula)
    return 0 if m is None else int(m.group(1) or 1)


def settings(sysdir):
    txt = (sysdir / "settings.yaml").read_text()
    z = int(re.search(r"charge carried:\s*(\d+)", txt).group(1))
    m = re.search(r"^\s+(Li|Na|K|Mg):\s*(-?\d+\.\d+)", txt, re.M)
    return z, (float(m.group(2)) if m else None)


def qe_energies(sysdir):
    """config name -> energy (eV, ASE constant) from relaxed_summary.csv"""
    out = {}
    for r in csv.DictReader((sysdir / "QE" / "relaxed_summary.csv").open()):
        if r["energy_Ry"]:
            out[r["config"]] = (float(r["energy_Ry"]) * RY, int(r["n_atoms"]), r["relax_finished"] == "True")
    return out


def states(sysdir, ion):
    """formula -> (n_ion, energy, source): lowest finished QE relaxation of that composition,
    plus VP/states entries that have no QE output at all (kept as they are)."""
    qe = {(r["state_file"], r["atoms_id"]): r["config"]
          for r in csv.DictReader((sysdir / "QE" / "states_crosscheck.csv").open())}
    E = qe_energies(sysdir)
    res, log = {}, []
    for st in sorted((sysdir / "VP" / "states").glob("*.csv")):
        formula = st.stem.split("-", 1)[1]
        cands = [(e, f"QE/{c}") for c, (e, n, fin) in E.items() if fin and c.rsplit("_", 1)[-1] == formula]
        for r in list(csv.reader(st.open()))[1:]:
            if len(r) >= 2 and qe.get((st.name, r[0].strip()), "") not in E:
                cands.append((float(r[1]), f"VP/states id {r[0].strip()} (no QE output)"))
        e, src = min(cands)
        if "VP/states" in src:
            log.append(f"{formula}: {src}")
        res[formula] = (n_ion(formula, ion), e, src)
    return res, log


def hull_and_voltage(st, mu, z):
    N = max(v[0] for v in st.values())
    full = next(v[1] for v in st.values() if v[0] == N)
    empty = next(v[1] for v in st.values() if v[0] == 0)
    pts = sorted((n / N, (e - (n / N) * full - (1 - n / N) * empty) / N, n, e, f, src)
                 for f, (n, e, src) in st.items())
    # lower convex hull (monotone chain)
    h = []
    for p in pts:
        while len(h) >= 2 and (h[-1][0] - h[-2][0]) * (p[1] - h[-2][1]) - (h[-1][1] - h[-2][1]) * (p[0] - h[-2][0]) <= 0:
            h.pop()
        h.append(p)
    volt = [(a[0], b[0], -(b[3] - a[3] - (b[2] - a[2]) * mu) / ((b[2] - a[2]) * z)) for a, b in zip(h, h[1:])]
    return pts, h, volt


def avg(v):
    return sum((b - a) * V for a, b, V in v) / sum(b - a for a, b, _ in v)


root = Path(sys.argv[1]).resolve()
out = root / "recomputed"
out.mkdir(exist_ok=True)
comp = [["system", "variant", "mu_eV", "x1", "x2", "V", "avg_V"]]
systems = sorted(root.glob("DFT-*"))
E_li = qe_energies(root / "DFT-Li2MnO3")
li_ref = E_li["Li"][0] / E_li["Li"][1]
for sysdir in systems:
    name = sysdir.name[4:]
    ion = ION_OF[name]
    z, mu_set = settings(sysdir)
    E = qe_energies(sysdir)
    mu_qe = E[ion][0] / E[ion][1] if ion in E and E[ion][2] else (li_ref if ion == "Li" else None)
    mu_set = mu_set if mu_set is not None else mu_qe  # Li2MnO3: settings point to Li.cif
    vp = [tuple(map(float, r)) for r in list(csv.reader((sysdir / "VP" / "voltage.csv").open()))[1:] if r]
    for a, b, V in vp:
        comp.append([name, "original", f"{mu_set:.4f}", f"{a:.4f}", f"{b:.4f}", f"{V:.4f}", f"{avg(vp):.4f}"])
    st, log = states(sysdir, ion)
    for v, mu in (("A_qe_energies", mu_set), ("B_qe_energies+qe_mu", mu_qe)):
        pts, h, volt = hull_and_voltage(st, mu, z)
        for a, b, V in volt:
            comp.append([name, v, f"{mu:.4f}", f"{a:.4f}", f"{b:.4f}", f"{V:.4f}", f"{avg(volt):.4f}"])
    # files for variant B (best estimate)
    with (out / f"{name}_voltage.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["x1", "x2", "V"])
        w.writerows([[f"{a:.6f}", f"{b:.6f}", f"{V:.6f}"] for a, b, V in volt])
    with (out / f"{name}_hull.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["formula", "n_ions", "x", "formation_energy_eV", "E_eV", "on_hull", "energy_source"])
        hs = {p[4] for p in h}
        w.writerows([[p[4], p[2], f"{p[0]:.6f}", f"{p[1]:.6f}", f"{p[3]:.6f}", p[4] in hs, p[5]] for p in pts])
    print(f"{name}: mu settings {mu_set:.4f}, mu QE {mu_qe:.4f}" + (" | kept from VP: " + "; ".join(log) if log else ""))

with (out / "voltage_comparison.csv").open("w", newline="") as f:
    csv.writer(f).writerows(comp)
