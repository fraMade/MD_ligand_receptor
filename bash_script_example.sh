#!/bin/bash

#SBATCH --job-name pdb_analysis
#SBATCH --time=24:00:00
#SBATCH -n 5
#SBATCH --account=account_X
#SBATCH --partition=g100_usr_prod
#SBATCH --error MD_lig_rcpt.err
#SBATCH --output MD_lig_rcpt.out

# Load the required modules for the CINECA G100 environment
# Attention to possible conflicts with others mpi libraries already imported
module load profile/lifesc
module load autoload gromacs/2021.2 
module load openmpi
module load python

# Insert you own virtual environment
source /g100_scratch/userexternal/mrX/my_venv/bin/activate


# The program will run 5 processes that will partition the first 100 picoseconds and analyze the ligand-protein interactions,
# the results will be save inside 'md_example/results'.

# Insert you own account and other details, e.g. the sbatch directives above
srun -A account_X -p g100_usr_prod -n 5 python MD-ligand-receptor.py -s ./md_example/molecular_dynamics.tpr -f ./md_example/molecular_dynamics.xtc -b 0 -e 100 -o ./md_example/ -l EHD

# Deactivate the virtual environment
deactivate