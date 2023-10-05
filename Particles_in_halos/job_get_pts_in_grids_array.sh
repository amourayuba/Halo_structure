#!/bin/bash
#SBATCH --account=def-taylor
#SBATCH -t 0-1:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=210G
#SBATCH --nodes=1
#SBATCH --array=0-4
source /home/ayuba/projects/def-taylor/ayuba/halo_profiles/bin/activate
python get_particles_in_grids.py $SLURM_ARRAY_TASK_ID


