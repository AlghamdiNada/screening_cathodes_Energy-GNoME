# DFT voltage profiles of existing cathodes (validation set)

0 K voltage profiles of six known cathodes, used to validate the computational screening method.
Each state along the (de)intercalation path is relaxed with Quantum ESPRESSO (`vc-relax`), the
convex hull is built, and voltage steps are taken between stable points.

This folder was cleaned on 2026-11-11 from the raw cluster working directories. Final inputs and
results are kept as they were (byte-identical), pw.x runs were collected per configuration, and
everything else was compressed into `archive/`, which is **not in git**; see
[Files not in git](#files-not-in-git).

## Method
- Code: Quantum ESPRESSO `pw.x` (GPU build, CINECA Leonardo), driven by
  `MC_MACE/mc_mace/pymc.py zerok-voltage settings.yaml` (see `job.sh`).
- Functional: r2SCAN via libxc, `input_dft = 'XC-000I-000I-000I-000I-497L-498L'`.
- Pseudopotentials: ONCV PBE-1.2 (`PP-r2scan` directory on the cluster).
- Relaxation: `vc-relax`, FIRE2, fmax = 0.05; etot/forc thresholds in each `settings.yaml`.
- Ion reference μ: energy per atom of the relaxed bulk metal (`QE/<Li|Na|K|Mg>/`), given in
  `settings.yaml` as a number, or as a CIF that the driver relaxes (Li2MnO3).

## Systems

| Folder | Ion | Full state | ecutwfc / ecutrho (Ry) | k-points | μ ion (eV/atom) | Steps | V range (V) | Avg V (V) | Final run |
|---|---|---|---|---|---|---|---|---|---|
| DFT-LiCoO2  | Li⁺  | Co6Li6O12   | 70 / 280 | 8×8×2 | −194.3176 ⚠ | 6 | 4.56 → 2.55 | 3.84 | 2025-05 |
| DFT-NaCoO2  | Na⁺  | Co6Na6O12   | 80 / 320 | 3×5×1 | −1161.234 | 4 | 3.82 → 1.47 | 2.60 | 2025-05 |
| DFT-LiFePO4 | Li⁺  | Fe8Li8O32P8 | 90 / 360 | 1×3×1 | −194.356  | 4 | 3.58 → 1.94 | 2.94 | 2025-07 |
| DFT-Li2MnO3 | Li⁺  | Li8Mn4O12   | 85 / 340 | 7×3×7 | from `Li.cif` (−194.357) | 6 | 4.99 → 2.01 | 3.70 ⚠ | 2025-07 |
| DFT-MgV2O4  | Mg²⁺ | Mg8O32V16   | 90 / 360 | 1×3×1 | −1481.803 | 5 | 1.00 → 0.03 | 0.64 | 2025-06 |
| DFT-KVPO4F  | K⁺   | F8K8O32P8V8 | 90 / 360 | 2×4×2 | −773.877  | 8 | 3.46 → 1.92 | 2.72 ⚠ | 2025-11 |

Avg V is the x-weighted mean of the steps in `VP/voltage.csv`. ⚠ = see
[Data-consistency findings](#data-consistency-findings).

## Layout

```
DFT-existing-cathodes-V-profile/
├── README.md
├── tools/                             scripts that built this folder (see below)
├── archive/DFT-<system>_old_runs.tar.gz
└── DFT-<system>/
    ├── *.cif                          structure (and ion reference) CIFs
    ├── settings.yaml                  driver input (relative paths: re-runnable from this folder)
    ├── full_simulation_settings.yaml  expanded settings written by the driver
    ├── job.sh                         SLURM script (Leonardo)
    ├── simulation.log                 log of the final run
    ├── VP/
    │   ├── voltage.csv                voltage steps: x1, x2, V
    │   ├── hull.csv                   formula, n ions, x, formation energy
    │   ├── thermo.csv, trj.xyz        run history
    │   └── states/                    NNN-<formula>.csv/.xyz: every configuration of each state
    ├── QE/                            one folder per relaxed configuration (see "QE folders")
    ├── plots/                         figures (KVPO4F, Li2MnO3)
    └── analysis/                      notebooks / helper scripts (KVPO4F, LiCoO2, NaCoO2)
```

`DFT-LiCoO2/analysis/pbe_settings.yaml` is an alternative PBE input, not the one used here.

## QE folders

`QE/<id>_<formula>/espresso.{pwi,pwo}` holds one pw.x `vc-relax` per configuration. `<id>` is the
`atoms_id` in `VP/states/*.csv`; `-1` / `[-1]` is a state with a single configuration (the full
state). `QE/<Li|Na|K|Mg>/` is the bulk-metal reference. Restarted runs keep the earlier outputs as
`N_espresso.{pwi,pwo}`. The relaxed cell is the **last `CELL_PARAMETERS` block** in the output,
followed by its `ATOMIC_POSITIONS`.

The driver overwrote `QE/` on every restart, so most runs were only saved inside the `QE_*.bak`
backups. `QE/` was rebuilt from those backups plus the last `QE/`: for each configuration, the most
recent run that finished (`JOB DONE`), or the latest attempt if none did. LiFePO4 `attempt_1/` is a
separate earlier attempt and was not used.

| System | Configs | Relaxation finished | Not finished |
|---|---|---|---|
| KVPO4F  | 39 | 39 | – |
| Li2MnO3 | 38 | 38 | – |
| LiCoO2  | 25 | 22 | `-1_Co27Li27O54` (separate 27-f.u. supercell test, never finished); `Li`, `backup-Li` (Li reference hit the max ionic steps) |
| LiFePO4 | 37 | 37 | – |
| MgV2O4  | 38 | 38 | – |
| NaCoO2  | 24 | 21 | `0_Co6Na5O12`, `[-1]_Co6Na6O12` (no complete step), `10_Co6NaO12` (stopped mid-relaxation) |

KVPO4F `1_F8K6O32P8V8` reached "End final coordinates" but has a CRASH file (bfgs "dE0s is
positive" in the final scf); its relaxed structure is usable.

Files added to each `QE/`:

| File | Content |
|---|---|
| `SOURCES.tsv` | which backup each config folder came from, and complete/INCOMPLETE |
| `<config>/relaxed.extxyz` | relaxed structure (last CELL_PARAMETERS + ATOMIC_POSITIONS); `info`: `energy_Ry`, `energy_eV`, `converged` (= "End final coordinates" reached), `source`, `config` |
| `relaxed_summary.csv` | one row per config: n_atoms, energy (Ry, eV), volume, relaxation finished, source pwo |
| `states_crosscheck.csv` | each `VP/states` entry compared with the QE energy of the same config |

Energies in eV use ASE's Ry = 13.6056930 eV. ASE's own `espresso-out` reader fails on these
outputs, so `tools/extract_relaxed.py` parses them directly.

Load the relaxed structures:
```python
from ase.io import read
atoms = read("DFT-LiCoO2/QE/0_Co6Li3O12/relaxed.extxyz")
atoms.info["energy_eV"], atoms.get_volume()
```

## Data-consistency findings

Of 191 state energies in `VP/states`, 185 match a QE output within 0.01 eV (using either constant
in point 1). The bulk-metal μ values in `settings.yaml` match the reference runs in `QE/`, except
for LiCoO2. `VP/` was **not** changed; corrected profiles are in [Recomputed voltages](#recomputed-voltages).

1. **Two Ry→eV constants are mixed in `VP/states`.** Most energies use 13.6056930 (ASE). 40 entries
   use 13.6057040: KVPO4F K7–K1 (33), MgV2O4 (6) and LiFePO4 Li5 id 28 (1). The KVPO4F and MgV2O4
   ones are configs that were re-run by hand (their folder has its own `job.sh`). That is
   0.03–0.05 eV per cell. In KVPO4F the K8 and K0 endpoints use the ASE value and the intermediate
   states use the other, so the first (K8→K7) and last (K1→K0) steps are off by about 0.03 V.
2. **Stale `VP/states` entries (the voltages are not affected).** In two cases `VP/states` holds an
   energy that the driver did *not* use. `VP/hull.csv` and `simulation.log` show that the hull and
   `voltage.csv` used the relaxed value:
   - LiFePO4 Fe8Li5O32P8 id 1: `states` has −42113.037 eV, 0.45 eV below any output. The hull used
     the surviving run (−42112.59 eV, the same as id 28).
   - Li2MnO3 LiMn4O12 (id −1): `states` has −16238.239 eV, an early unrelaxed step of
     `0_LiMn4O12/1_espresso.pwo`. The hull used the relaxed value (0.36 eV lower).
3. **LiCoO2 μ_Li = −194.3176 eV** doesn't come from any surviving Li run. The two LiCoO2 Li runs
   give −194.086 and −194.343 (neither finished). LiFePO4 and Li2MnO3 use −194.356/−194.357. With
   that value, all LiCoO2 voltages are 0.039 V lower (see [Recomputed](#recomputed-voltages)).
4. **NaCoO2 Na6 (full state) and Na5 id 0**: the energies in `VP` have no surviving QE output.
5. Minor: KVPO4F K7 id 1 uses an intermediate step of `1_espresso.pwo` (0.012 eV from the final
   value). Li2MnO3 Mn4O12 uses the single-scf restart `espresso.pwo`, while `relaxed.extxyz` comes
   from the full relaxation in `1_espresso.pwo`.

## Recomputed voltages

`tools/recompute_voltage.py` reimplements the driver's hull and voltage step: formation energy per
host relative to the full and empty states, lower convex hull, and
V = −[E(n₂) − E(n₁) − (n₂−n₁)·μ] / ((n₂−n₁)·z). Fed the `VP/states` energies, it reproduces
`VP/voltage.csv` for KVPO4F, LiCoO2, MgV2O4 and NaCoO2. Li2MnO3 and LiFePO4 differ only because of
the stale entries in point 2. It writes to `recomputed/`:

- **A**: each state's energy is the lowest *finished* QE relaxation of that composition, all
  converted with one constant (ASE). μ is as in `settings.yaml`. `VP/states` entries with no QE
  output are kept; this applies only to NaCoO2 Na6.
- **B** (best estimate): A, plus μ from a finished bulk-metal QE run with the same constant.
  LiCoO2 and LiFePO4 have no finished Li run of their own, so they use Li2MnO3's Li run
  (−194.3567 eV). The Li runs use slightly different cutoffs (70–90 Ry); the size of that effect
  hasn't been checked.

| System | Avg V original | Avg V A | Avg V B | Max step change A | Max step change B | Same hull points |
|---|---|---|---|---|---|---|
| LiCoO2  | 3.836 | 3.836 | **3.797** | 0.000 | 0.039 | yes |
| NaCoO2  | 2.598 | 2.598 | 2.598 | 0.003 | 0.003 | yes |
| LiFePO4 | 2.941 | 2.941 | 2.940 | 0.000 | 0.001 | yes |
| Li2MnO3 | 3.697 | 3.696 | 3.696 | 0.024 | 0.024 | yes |
| MgV2O4  | 0.641 | 0.641 | 0.641 | 0.023 | 0.023 | yes |
| KVPO4F  | 2.717 | 2.717 | 2.717 | 0.033 | 0.032 | yes |

The hull (which states are stable) is unchanged everywhere. The individual steps that move are:
- KVPO4F: the first step drops 3.457 → 3.425 V; the steps at x = 0.75–0.875 and 0.875–1 rise
  2.077 → 2.096 and 1.923 → 1.940 V.
- MgV2O4: the middle steps shift by ±0.02 V.
- Li2MnO3: the first two steps shift 4.989/4.110 → 4.965/4.125 V.
- LiCoO2: every step drops by 0.039 V (μ only).

The mixed conversion constants (point 1) only move individual steps and leave the average
voltages unchanged; μ_Li moves the LiCoO2 average.

Files: `recomputed/<system>_voltage.csv` and `<system>_hull.csv` (variant B; the hull lists the
energy source of each state), and `recomputed/voltage_comparison.csv` (all variants, all steps).

## Files not in git

The pw.x outputs (`QE/**/espresso.pwo` and restart `N_espresso.pwo`: 230 files, 280 MB raw,
largest 19 MB) **are** in git. `.gitignore` excludes:
- `archive/`: the old runs and backups (200 MB, see below; one file is over GitHub's 100 MB
  limit). Kept locally and on the cluster only.
- SLURM `*.out/*.err`, `espresso.err`, `CRASH` inside `QE/`.
- `*.bak*` backups and OS/editor junk.

## Archive
Each `archive/DFT-<system>_old_runs.tar.gz` (local only) holds, with original paths:
the automatic `*.bak` restart backups (QE, plots, logs, VP files), `convergence/` tests
(ecutwfc / k-point), `lowest_energy_states/`, earlier attempts (`LiFePO4/attempt_1`, nested
`Li2MnO3/Li2MnO3`), manual state copies (`states-b`, `states (Copy)`, `my-back-states*`) and
old SLURM `.out/.err` files. The runs collected into `QE/` are also still in the archive.
macOS/vim junk (`.DS_Store`, `._*`, `*.swp`) was dropped.

Extract one system:
```bash
tar -xzf archive/DFT-LiCoO2_old_runs.tar.gz            # -> DFT-LiCoO2/... (no file overlaps)
tar -tzf archive/DFT-LiCoO2_old_runs.tar.gz | less     # list only
```

## tools/

| Script | Purpose |
|---|---|
| `make_clean.py <orig> <clean>` | sort files into kept / `analysis/` / archive / dropped; write the archives |
| `collect_qe.py <orig sys> <clean sys>` | rebuild `QE/` from `QE_*.bak` + `QE/`, write `SOURCES.tsv` |
| `extract_relaxed.py <clean>` | write `relaxed.extxyz`, `relaxed_summary.csv`, `states_crosscheck.csv` |
| `recompute_voltage.py <clean>` | recompute hull + voltages (variants A/B), write `recomputed/` |
| `plot_voltage_comparison.py <clean>` | plot original vs recomputed profiles → `recomputed/voltage_comparison.{png,pdf}` |

`make_clean.py` and `collect_qe.py` were run once on the raw cluster folder (no longer present) and
are kept as a record. `extract_relaxed.py` needs the `.pwo` files. The last two scripts run on
what is in git.
