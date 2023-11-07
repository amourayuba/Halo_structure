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
Given the particle files produced above and located in "folder", and AHF halos files located in folder/AHF/halos/* 
this will associate particles and halos. 

Things to potentially modify in the code: 
sims: list of the simulation names. 
snp: snapshot number to consider 
folder: where are the gadget files located, and where the output files will be stored
divisions: grid size. Number of cells = divisions**3. 

All of the above have to be coherent with the grid files generated from get_particles_in_grids.py. 
nlim how many virial radii away from the halo centre to gather particles

Will save a list where each element has halo_id, and particle information (id, velocity, position) for each particle
in the halo

#### *) job_pts_in_halos_array.sh
sbatch script to run on compute canada. Array job, parameter is the list of simulation indices. 
Can run with #SBATCH --array=x, y, z if one wants to do it for only simulation indices x, y, z
Second parameter is the minimum halo mass to consider. 

### 2) Measuring structural properties for each halos ./structural_properties

#### *) fits2d.py 
Functions calculating different halo properties, 3D or 2D. These will be used to 

#### *) get_2d_3d_props.py 
For a given simulation id, and mass lim, will calculate and save a set of 2D and 3D proprties of halos. 
First set the simulation folder in the variable folder

Files required in the folder "folder": 
{sim}_prefixes.txt : list of AHF file prefiexes
/folder/AHF/halos/*_AHF_halos : AHF halo files of the form {prefix}_AHF_halos
folder/parts_near_halos*.npy : particle in halos files generated through particles_in_halos scripts

outputs: 
###### props3D*.npy : 
list of 3D halo property, each element (represents a halo) has the following components: 
*) conc3d: concentration parameter found through the 3D fit of the density profile 
*) chirho3d: array of two elements. 1) result of fits2d chi_square function, 2 result of the log_chi_square function. Both 
for the 3d density profile, fitted vs real one.
*) chim3d: same as chirho3d for the mass profile instead of density profiles

###### props2D*.npy
list of 2D projected halo property, each element (represents a halo) has the following components:
*) conc2d: concentration parameter found through the 2D fit of the projected density profile 
*) chirho2d: array of two elements. 1) result of fits2d chi_square function, 2 result of the log_chi_square function. Both 
for the 2d projected density profile, fitted vs real one.
*) chim2d same as chirho2d for the mass profile instead of density profiles
*) axs2d 2d axis ratio between the second and first major axis, gives ellipticity in the x-y plane
*) vecs2d vectors major axes direction in the x-y plane
*) mboff2d distance between the most bound particle and density peak of the halo in the x-y plane
*) com2d distance between the centre of mass of the halo and peak of the density profile in the x-y plane

#### *) job_2d_3d_props.sh 
sbatch script to run on compute canada. Array job, parameter is the list of simulation indices. 
Can run with #SBATCH --array=x, y, z if one wants to do it for only simulation indices x, y, z
Second parameter is the minimum halo mass to consider. 
Need to set up a python environment, and install the relevent packages.

