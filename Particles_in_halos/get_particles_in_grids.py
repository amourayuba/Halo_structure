import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys, os
import readfof
import readgadget
import readsnap
import h5py
import csv
from scipy.optimize import minimize
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
from fits2ds import *


#################### FUNCTIONS #####################

# create 3D array, each element is a cell of the grid used to divide 
# the volume of sim. Each element has the coordinates of the particles 
# in them, as ordered in the pos array                    

def get_part_inGrid(divisions, BoxSize, pos):
    '''
    :param divisions: int, number of divisions for the grid. Ncells = division**3
    :param BoxSize: float, size of the cosmological box where particles are
    :param pos: array[N,3] 3D positions of the particles
    :return: numpy array [3*division, [Nparticles in division]], dtype = object. Particle indices that
    are in each cell.
    '''
    print('\nsorting particle positions')
    index = (pos * divisions / BoxSize).astype(
        np.uint32)  # This attributes a triplet between 0 and division to each particle.
    index[np.where(
        index == divisions)] = divisions - 1  # I'm not sure why, but it puts the particles on the edge, one grid before
    print(np.min(index), '< index <', np.max(index))

    cell = index[:, 0] * divisions ** 2 + index[:, 1] * divisions + index[:, 2]  # seems like a way to incode the
    # 3D index in 1D cell. It works because 0<index<divisions, so this is a bijective mapping
    del index
    print(np.min(cell), '< cell <', np.max(cell))

    sorted_ids = []
    for i in range(divisions ** 3):
        sorted_ids.append([])

    for i in range(len(cell)):  # len(cell) is the number of particles we have
        sorted_ids[cell[i]].append(i)  # for each cell, we put the id's of all particles in that cell, this only works
        # because the id's were sorted, so that i is the id of the particle.
    sorted_ids = np.array(sorted_ids)

    return sorted_ids


if __name__ == "__main__":
    sims = ['m4s7', 'm3s85', 'm25s9', 'm35s9', 'm35s7', 'm25s85', 'm2s8', 'm4s8', 'm2s9', 'm3s8_50', 'm3s8', 'm35s75',
            'm4s9', 'm3s9', 'm25s75', 'm2s1', 'm3s7', 'm3s75', 'm2s7', 'm25s8', 'm35s8', 'm3s8b', 'm35s85']

    simid = int(sys.argv[1])
    sim = sims[simid]  # 'm2s9'
    snp = 118
    folder = '/home/ayuba/scratch/{}/'.format(sim)
    print(folder)

    snap = 'snapdir_{:03}/snapshot_{:03}'.format(snp, snp)
    head = readgadget.header(folder + snap)
    BoxSize = head.boxsize  # Mpc/h
    Nall = head.nall  # Total number of particles
    Masses = head.massarr * 1e10  # Masses of the particles in Msun/h
    mpart = Masses[1]
    # compute mean density from snapshot
    mean_background = Nall[1] / BoxSize ** 3
    print(mean_background)

    # order the particle positions as the ID (position of particle with ID 1 is written in position 0 of this array, and so on)
    ptype = [1]
    save = True
    # sort the particle positions in cells
    bins = 50
    divisions = 200
    if save:
        print('Reading pt positions and ids')

        pos_0 = readgadget.read_block(folder + snap, "POS ", ptype).astype(np.float64)  # Mpc/h
        vel_0 = readgadget.read_block(folder + snap, "VEL ", ptype).astype(np.float64)  # km/s

        grid_0 = get_part_inGrid(divisions, BoxSize, pos_0)
        np.save(folder + '1024_{}_grid_0_div{}_snap{}.npy'.format(sim, divisions, snp), grid_0)
        np.save(folder + '1024_{}_pos_sorted_div{}_snap{}.npy'.format(sim, divisions, snp), pos_0)
        np.save(folder + '1024_{}_vel_sorted_div{}_snap{}.npy'.format(sim, divisions, snp), vel_0)

    else:
        print('Loading pt positions and ids')
        grid_0 = np.load(folder + '1024_{}_grid_0_div{}_snap{}.npy'.format(sim, divisions, snp), allow_pickle=True)
        pos_0 = np.load(folder + '1024_{}_pos_sorted_div{}_snap{}.npy'.format(sim, divisions, snp), allow_pickle=True)
        vel_0 = np.load(folder + '1024_{}_vel_sorted_div{}_snap{}.npy'.format(sim, divisions, snp), allow_pickle=True)
