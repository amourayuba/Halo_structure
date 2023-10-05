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
The set of codes in this folder allows to associate particles' positions and velocities to a set of halos. 

#### REQUIRED : 



