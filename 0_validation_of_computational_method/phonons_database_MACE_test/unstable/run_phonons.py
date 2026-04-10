#!/usr/bin/env python3
import os
import shutil

cwd = './'
os.chdir(cwd)

for files in os.listdir(f'{cwd}/{folder_name}'):
    if os.path.isdir(f'./{files}'):
        print(files)
        path = f'./{files}'
        code = f'''#!/usr/bin/env python3
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

info = [f'{files}']
calculator = MACECalculator(model_path='mc_mace/GNoME-screening/uMLFF_are_ready_for_phonons_data_MACE_test/mace-omat-0-medium.model', device='cuda')

if not os.path.exists(f'{cwd}/{path}/relaxed.cif'):
    atoms = read(f'{cwd}/{path}/{files}')
    atoms.calc = calculator
    print('Starting relaxation')
    dyn = FIRE2(atoms)
    dyn.run(fmax=1e-4)
    write(f'{cwd}/{path}/relaxed.cif', atoms)
else:
    atoms = read(f'{cwd}/{path}/relaxed.cif')

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
        ph.clean()
        print('Starting phonon calculation with 4x4x4 supercell')
        ph = Phonons(atoms, calculator, supercell=(4, 4, 4))
        ph.run()
        finished = True
        super_cell_dim = 4
    except Exception as e: print(e)
if not finished:
    try:
        ph.clean()
        print('Starting phonon calculation with 3x3x3 supercell')
        ph = Phonons(atoms, calculator, supercell=(3, 3, 3))
        ph.clean()
        ph.run()
        finished = True
        super_cell_dim = 3
    except Exception as e: print(e)
if not finished:
    try:
        ph.clean()
        print('Starting phonon calculation with 2x2x2 supercell')
        ph = Phonons(atoms, calculator, supercell=(2, 2, 2))
        ph.run()
        finished = True
        super_cell_dim = 2
    except Exception as e: print(e)
try:
    ph.read(method='frederiksen', acoustic=True)
    ph.clean()
    path_ = atoms.cell.bandpath()
    bs, modes = ph.get_band_structure(path_, modes=True)
    with open(f'{cwd}/{path}/phonon_bandstructure.json', 'w') as f:
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

    plt.savefig(f'{cwd}/{path}/bands.pdf')
    print('Phonon band structure calculation completed successfully')

except Exception as e: print(e)

# ---- load bandstructure (ASE JSON) ----
if os.path.exists(f'{cwd}/{path}/phonon_bandstructure.json'):
    with open(f'{cwd}/{path}/phonon_bandstructure.json', 'r') as f:
        bs = BandStructure.read(f) 
    n_negative_bands = len(np.where(bs.energies < -0.004)[0])
    exists_imaginary = n_negative_bands > 0
else:
    print(f'No phonon_bandstructure.json found in {cwd}/{path}/phonon_bandstructure.json, skipping...')

info.append(potentialenergy)
info.append(super_cell_dim)
info.append(atoms.get_cell().lengths()*super_cell_dim)
info.append(atoms.get_number_of_atoms()*super_cell_dim**3)
info.append(n_negative_bands)
info.append(exists_imaginary)

df = pd.DataFrame([info], columns=['id', 'chemical_formula', 'n_atoms', 'cell_lengths', 'potential_energy (eV)', 'supercell_dim', 'supercell_lengths', 'n_atoms_supercell', 'n_imaginary_bands', 'has_imaginary'])
df.to_csv(f'{cwd}/{path}/info.csv', index=False)

'''
        job = f'''#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --qos normal
#SBATCH --time=02:00:00
#SBATCH --gres=gpu:1
#SBATCH --job-name={files}       # nome della simulazione
##SBATCH --output=%x_%j.out
##SBATCH --error=%x_%j.err
#SBATCH -A xxxx
#SBATCH -p boost_usr_prod

module purge
module load programs/nvhpc-24.1/modulefiles/nvhpc-hpcx-cuda12/24.1
source volta/bin/activate

python3 run_phonons.py 
'''
        with open(f'{path}/run_phonons.py', 'w') as f:
            f.write(code)
        with open(f'{path}/job.sh', 'w') as f:
            f.write(job)
        print(f'Created folder and script for {files}')
