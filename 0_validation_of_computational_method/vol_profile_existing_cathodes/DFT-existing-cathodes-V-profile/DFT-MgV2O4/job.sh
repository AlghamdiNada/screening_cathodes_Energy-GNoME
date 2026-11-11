#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --qos normal         # boost_qos_dbg, normal
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:1
#SBATCH --job-name=MgV2O4        # nome della simulazione
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH -A INR24_CHIROC_0 #IscrB_NEXT-LIB
#SBATCH -p boost_usr_prod

module purge
module load /leonardo/pub/userexternal/nalghamd/programs/nvhpc-24.1/modulefiles/nvhpc-hpcx-cuda12/24.1
source /leonardo/pub/userexternal/nalghamd/env/volta/bin/activate
#module load profile/chem-phys
#module load quantum-espresso/7.2--openmpi--4.1.4--nvhpc--23.1-openblas-cuda-11.8

export PYTHONUNBUFFERED=TRUE
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#export CRAY_CUDA_MPS=1
#export ASE_ESPRESSO_COMMAND_BASE='srun --cpu-bind=cores pw.x'

#ulimit -s unlimited

#srun --cpu-bind=cores pw.x -in pwscf.in
#python3 leonardo_scratch/large/userexternal/nalghamd/mc_mace/MC_MACE/mc_mace/pymc.py zerok-voltage settings.yaml --log-level debug
python3 ../MC_MACE/mc_mace/pymc.py zerok-voltage settings.yaml --log-level debug
