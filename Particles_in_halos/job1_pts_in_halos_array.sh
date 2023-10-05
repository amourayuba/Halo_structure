#!/bin/bash
#SBATCH -t 0-1:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=310G
#SBATCH --nodes=1
#SBATCH --array=4,3,10

source /home/ayuba/projects/def-taylor/ayuba/halo_profiles/bin/activate

#Assuming the files for particles in grid exist, and the files in halos do too. first arg is simulation id, second is mlim 
python get_particles_in_halos.py $SLURM_ARRAY_TASK_ID 0


