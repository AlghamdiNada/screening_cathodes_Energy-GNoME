#!/usr/bin/env python3
from ase.io import read, write
from mace.calculators import MACECalculator
from ase.phonons import Phonons

from ase.optimize import BFGS, FIRE2
import numpy as np
from ase.spectrum.band_structure import BandStructure
import matplotlib.pyplot as plt
import os
import shutil
import pandas as pd

info = [f'c1551f877c']
calculator = MACECalculator(model_path='/leonardo_scratch/large/userexternal/nalghamd/mc_mace/GNoME-screening/correct_phonon_run/mace-omat-0-medium.model', device='cuda')

if not os.path.exists(f'/leonardo_scratch/large/userexternal/nalghamd/mc_mace/GNoME-screening/correct_phonon_run/Mg/c1551f877c/relaxed.cif'):
    atoms = read(f'/leonardo_scratch/large/userexternal/nalghamd/mc_mace/GNoME-screening/correct_phonon_run/Mg/c1551f877c/c1551f877c.CIF')
    atoms.calc = calculator
    print('Starting relaxation')
    dyn = FIRE2(atoms)
    dyn.run(fmax=1e-4)
    write(f'/leonardo_scratch/large/userexternal/nalghamd/mc_mace/GNoME-screening/correct_phonon_run/Mg/c1551f877c/relaxed.cif', atoms)
else:
    atoms = read(f'/leonardo_scratch/large/userexternal/nalghamd/mc_mace/GNoME-screening/correct_phonon_run/Mg/c1551f877c/relaxed.cif')

info.append(atoms.get_chemical_formula())
info.append(atoms.get_number_of_atoms())
info.append(atoms.get_cell().lengths())

atoms.calc = calculator

potentialenergy = atoms.get_potential_energy()
#ph.clean()

finished = False
super_cell_dim = 0

try:
    print('Starting phonon calculation with 5x5x5 supercell')
    ph = Phonons(atoms, calculator, supercell=(5, 5, 5))
    ph.run()
    finished = True
    super_cell_dim = 5
except Exception as e: print(e)
if not finished:
    try:
        print('Starting phonon calculation with 4x4x4 supercell')
        ph = Phonons(atoms, calculator, supercell=(4, 4, 4))
        ph.clean()
        ph.run()
        finished = True
        super_cell_dim = 4
    except Exception as e: print(e)
if not finished:
    try:
        print('Starting phonon calculation with 3x3x3 supercell')
        ph = Phonons(atoms, calculator, supercell=(3, 3, 3))
        ph.clean()
        ph.run()
        finished = True
        super_cell_dim = 3
    except Exception as e: print(e)
if not finished:
    try:
        print('Starting phonon calculation with 2x2x2 supercell')
        ph = Phonons(atoms, calculator, supercell=(2, 2, 2))
        ph.clean()
        ph.run()
        finished = True
        super_cell_dim = 2
    except Exception as e: print(e)
try:
    ph.read(method='frederiksen', acoustic=True)
    ph.clean()
    path_ = atoms.cell.bandpath()
    bs, modes = ph.get_band_structure(path_, modes=True)
    with open(f'/leonardo_scratch/large/userexternal/nalghamd/mc_mace/GNoME-screening/correct_phonon_run/Mg/c1551f877c/phonon_bandstructure.json', 'w') as f:
        bs.write(f)
    fc = ph.get_force_constant()  # real-space force constants array
    np.save('force_constants.npy', fc)
    dos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=1000, width=1e-6)
    energies = dos.get_energies()
    weights  = dos.get_weights()
    # Save DOS as a compact NumPy .npz (easy to load anywhere)
    np.savez('phonon_dos.npz', energies=energies, dos=weights, info=dos.info)

    fig = plt.figure(figsize=(7, 4))
    ax = fig.add_axes([0.12, 0.07, 0.67, 0.85])
    emax = 0.035
    bs.plot(ax=ax, emin=-emax, emax=emax)

    dosax = fig.add_axes([0.8, 0.07, 0.17, 0.85])
    dosax.fill_between(
        dos.get_weights(),
        dos.get_energies(),
        y2=0,
        color='grey',
        edgecolor='k',
        lw=1,
    )
    dosax.set_ylim(-emax, emax)
    dosax.set_yticks([])
    dosax.set_xticks([])
    dosax.set_xlabel('DOS')

    plt.savefig(f'/leonardo_scratch/large/userexternal/nalghamd/mc_mace/GNoME-screening/correct_phonon_run/Mg/c1551f877c/bands.pdf')
    print('Phonon band structure calculation completed successfully')

except Exception as e: print(e)

# ---- load bandstructure (ASE JSON) ----
if os.path.exists(f'/leonardo_scratch/large/userexternal/nalghamd/mc_mace/GNoME-screening/correct_phonon_run/Mg/c1551f877c/phonon_bandstructure.json'):
    with open(f'/leonardo_scratch/large/userexternal/nalghamd/mc_mace/GNoME-screening/correct_phonon_run/Mg/c1551f877c/phonon_bandstructure.json', 'r') as f:
        bs = BandStructure.read(f) 
    n_negative_bands = len(np.where(bs.energies < -0.004)[0])
    exists_imaginary = n_negative_bands > 0
else:
    print(f'No phonon_bandstructure.json found in /leonardo_scratch/large/userexternal/nalghamd/mc_mace/GNoME-screening/correct_phonon_run/Mg/c1551f877c/phonon_bandstructure.json, skipping...')

info.append(potentialenergy)
info.append(super_cell_dim)
info.append(atoms.get_cell().lengths()*super_cell_dim)
info.append(atoms.get_number_of_atoms()*super_cell_dim**3)
info.append(n_negative_bands)
info.append(exists_imaginary)

df = pd.DataFrame([info], columns=['id', 'chemical_formula', 'n_atoms', 'cell_lengths', 'potential_energy (eV)', 'supercell_dim', 'supercell_lengths', 'n_atoms_supercell', 'n_imaginary_bands', 'has_imaginary'])
df.to_csv(f'/leonardo_scratch/large/userexternal/nalghamd/mc_mace/GNoME-screening/correct_phonon_run/Mg/c1551f877c/info.csv', index=False)

