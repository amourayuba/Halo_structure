# Halo_structure

## Gets structural properties of halos from gadget simulation outputs (particles) and Amiga Halo Finder outputs (halos). This is achieved through several steps. 
All of the codes below assume that files and subfolders will be located in a directory "folder" (name of the python variable) which will have the structure '*/{sim}/' where sim is going to be 
the name of the simulation considered. Feel free to modify the python variable in "folder" in each of the python codes below as needed. All files produces with the codes below are stored by default in 
"folder". 
### IMPORTANT: By default
AHF files need to be present in folder/AHF/*AHF_halos
Gadget snapshots have to be present in folder/snapdir_*/*.hdf5 
In order to read AHF files, a {sim}_prefixes.txt doc is also located in folder, each line of this text file is the AHF prefix output for each snapshot, starting from the most recent.
i.e. each AHF file will have the form prefix.AHF_{halos, particles, profiles, substructure} 

### 1) Associating particles to halos ./Particles_in_halos 
The set of codes in this folder allows us to associate particles' positions and velocities to a set of halos. 

#### REQUIRED : 
numpy, Pylians : https://pylians3.readthedocs.io/en/master/, pandas

As well as the requirements above for gadget and AHF files to be located in "folder". 

#### *) get_particles_in_grids.py 
The first step is to put particles in a grid of chosen size "divisions". 
Things to potentially modify in the code: 
sims: list of the simulation names. 
snp: snapshot number to consider 
folder: where are the gadget files located, and where the output files will be stored
divisions: grid size. Number of cells = divisions**3. 
This will produce 3 files where each is a list of particle_id, particle_positions, particle_velocities in each cell. Index of a cell is the index of the list.

#### *) job_get_pts_in_grids_array.sh 
sbatch script to run on compute canada. Array job, parameter is the list of simulation indices. Can run with #SBATCH --array=x, y, z if one wants to do it for only simulation indices x, y, z

#### *) get_particles_in_halos.py

Given the particle files produced above
