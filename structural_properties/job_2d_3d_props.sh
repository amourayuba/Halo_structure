#!/bin/bash
#SBATCH -t 1-20:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=310G
#SBATCH --nodes=1
#SBATCH --array=0-19

source /home/ayuba/projects/def-taylor/ayuba/halo_profiles/bin/activate
python get_2d_3d_props.py $SLURM_ARRAY_TASK_ID 0


