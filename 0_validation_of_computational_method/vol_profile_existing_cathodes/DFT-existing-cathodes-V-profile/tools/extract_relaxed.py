"""For each clean DFT-<sys>/QE/<config>/, write relaxed.extxyz (last frame of the pw.x output)
and QE/relaxed_summary.csv; then cross-check energies against VP/states/*.csv."""
import sys, csv, re
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import write
from ase.units import Bohr

RY_EV = 13.605693012183622  # ASE CODATA-2014 Ry->eV (ase.units.Ry)
RY_ALT = 13.605703976       # older constant found in some VP/ energies

root = Path(sys.argv[1]).resolve()


def parse_pwo(pwo):
    """Last CELL_PARAMETERS + following ATOMIC_POSITIONS block and last '!' total energy (eV)."""
    lines = pwo.read_text(errors="ignore").splitlines()
    alat = None
    for l in lines:
        if "lattice parameter (alat)" in l:
            alat = float(l.split("=")[1].split()[0]) * Bohr  # Angstrom
            break
    ic = max(i for i, l in enumerate(lines) if l.startswith("CELL_PARAMETERS"))
    hdr = lines[ic]
    cell = np.array([[float(x) for x in lines[ic + k].split()[:3]] for k in (1, 2, 3)])
    if "angstrom" in hdr:
        pass
    elif "bohr" in hdr:
        cell *= Bohr
    else:  # (alat= X) in bohr
        m = re.search(r"alat\s*=\s*([\d.]+)", hdr)
        cell *= float(m.group(1)) * Bohr if m else alat
    ip = next(i for i in range(ic, len(lines)) if lines[i].startswith("ATOMIC_POSITIONS"))
    phdr = lines[ip]
    syms, pos = [], []
    for l in lines[ip + 1:]:
        p = l.split()
        if len(p) < 4 or not re.match(r"^[A-Za-z]", p[0]):
            break
        try:
            xyz = [float(x) for x in p[1:4]]
        except ValueError:
            break
        syms.append(re.sub(r"\d+$", "", p[0]))
        pos.append(xyz)
    pos = np.array(pos)
    if "crystal" in phdr:
        a = Atoms(syms, scaled_positions=pos, cell=cell, pbc=True)
    else:
        if "bohr" in phdr:
            pos *= Bohr
        elif "alat" in phdr:
            pos *= alat
        a = Atoms(syms, positions=pos, cell=cell, pbc=True)
    en = [l for l in lines if l.startswith("!") and "total energy" in l]
    if not en:
        raise ValueError("no converged scf energy")
    return a, float(en[-1].split("=")[1].split()[0])  # Ry


def pick_pwo(c):
    pwos = list(c.glob("*espresso.pwo"))
    if not pwos:
        return None, False, "no pwo"
    fin = {p: "End final coordinates" in p.read_text(errors="ignore") for p in pwos}
    main = c / "espresso.pwo"
    if fin.get(main):
        return main, True, ""
    num = sorted((p for p in pwos if fin[p] and p.name[0].isdigit()),
                 key=lambda p: int(p.name.split("_")[0]))
    if num:
        return num[-1], True, ""
    return (main if main.exists() else sorted(pwos)[-1]), False, "relaxation not finished"


mismatches = []
for sysdir in sorted(root.glob("DFT-*")):
    qe = sysdir / "QE"
    status = {}
    src = qe / "SOURCES.tsv"
    if src.exists():
        for row in list(csv.reader(src.open(), delimiter="\t"))[1:]:
            status[row[0]] = row[1]
    rows = []
    for c in sorted(p for p in qe.iterdir() if p.is_dir()):
        pwo, conv, note = pick_pwo(c)
        base = [c.name, c.name.rsplit("_", 1)[-1]]
        if pwo is None:
            rows.append(base + ["", "", "", "", False, status.get(c.name, ""), "", note]); continue
        try:
            a, e_ry = parse_pwo(pwo)
        except Exception as ex:
            rows.append(base + ["", "", "", "", False, status.get(c.name, ""), pwo.name,
                                "no complete ionic step in output"]); continue
        e = e_ry * RY_EV
        a.info.update(energy_eV=e, energy_Ry=e_ry, converged=conv, source=pwo.name, config=c.name)
        write(c / "relaxed.extxyz", a, format="extxyz")
        rows.append(base + [len(a), f"{e_ry:.8f}", f"{e:.6f}", f"{a.get_volume():.4f}",
                            conv, status.get(c.name, ""), pwo.name, note])
    with (qe / "relaxed_summary.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["config", "formula", "n_atoms", "energy_Ry", "energy_eV", "volume_A3",
                    "relax_finished", "run_status", "source_pwo", "note"])
        w.writerows(rows)
    by_name = {r[0]: r for r in rows if r[4]}

    # cross-check against VP/states (energies there are printed to ~0.01 eV)
    TOL = 0.01
    out = [["state_file", "atoms_id", "E_states_eV", "config", "E_QE_eV", "diff_eV", "result",
            "closest_config_same_formula", "closest_diff_eV"]]
    counts = {"match": 0, "match (Ry=13.605704)": 0, "MISMATCH": 0, "NO_OUTPUT": 0}
    for st in sorted((sysdir / "VP" / "states").glob("*.csv")):
        formula = st.stem.split("-", 1)[1]
        same = [r for r in by_name.values() if r[1] == formula]
        for r in [x for x in list(csv.reader(st.open()))[1:] if len(x) >= 2]:
            sid, e_ref = r[0].strip(), float(r[1])
            cands = [f"{sid}_{formula}", f"[{sid}]_{formula}"]
            if sid == "-1" and len(same) == 1:  # full/empty state saved under another name
                cands.append(same[0][0])
            hit = next((by_name[n] for n in cands if n in by_name), None)
            close = min(same, key=lambda q: abs(float(q[4]) - e_ref), default=None)
            cl = [close[0], f"{float(close[4]) - e_ref:+.4f}"] if close else ["", ""]
            if hit is None:
                res, row = "NO_OUTPUT", ["", "", ""]
            else:
                e_ry = float(hit[3])
                d = e_ry * RY_EV - e_ref
                if abs(d) <= TOL:
                    res = "match"
                elif abs(e_ry * RY_ALT - e_ref) <= TOL:
                    res = "match (Ry=13.605704)"
                else:
                    res = "MISMATCH"
                row = [hit[0], hit[4], f"{d:+.4f}"]
            counts[res] += 1
            out.append([st.name, sid, f"{e_ref:.4f}"] + row + [res] + cl)
            if not res.startswith("match"):
                mismatches.append(f"{sysdir.name}: {st.name} id {sid}: {res} {row[0]} {row[2]} | closest {cl[0]} {cl[1]}")
    with (qe / "states_crosscheck.csv").open("w", newline="") as f:
        csv.writer(f).writerows(out)
    nconv = sum(1 for r in rows if r[6] is True)
    print(f"{sysdir.name}: {len(rows)} configs, {nconv} relaxation finished, states check {counts}")
print("\n".join(mismatches) or "no mismatches")
