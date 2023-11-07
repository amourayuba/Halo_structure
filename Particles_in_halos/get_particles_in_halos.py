import numpy as np
import sys
import pandas as pd
import readgadget


def get_part_nearHalo(pos_h, dist_max, divisions, BoxSize, sorted_ids):
    '''
    :param pos_h: array[3,1], float center of halo
    :param dist_max: float. maximum distance around the halo centre up to which to consider particles
    :param divisions: int. number of divisions for the grid. Ncells = division**3
    :param BoxSize: float, size of the cosmological box where particles are
    :param sorted_ids: numpy array [3*division, [Nparticles in division]], dtype = object. Particle indices that.
    Output of get_part_inGrid()
    :return: numpy array[N,1] int. Indices of all particles within dist_max of pos_h
    '''

    # which cells I need to take into account
    index_min = np.floor(((pos_h - dist_max) * divisions / BoxSize))
    index_max = np.floor(((pos_h + dist_max) * divisions / BoxSize))
    index_min = index_min.astype(np.int32)
    index_max = index_max.astype(np.int32)

    i_num = np.arange(index_min[0], index_max[0] + 1) % divisions
    j_num = np.arange(index_min[1], index_max[1] + 1) % divisions
    k_num = np.arange(index_min[2], index_max[2] + 1) % divisions

    # count the particles in the selected cells                        
    length = 0
    par_indexes = []
    for i in range(len(i_num)):
        for j in range(len(j_num)):
            for k in range(len(k_num)):  # go through all cells, count nparticles in each cell
                box_id = divisions ** 2 * i_num[i] + divisions * j_num[j] + k_num[k]
                length += len(sorted_ids[box_id])
    par_indexes = np.empty(length, dtype=np.uint32)

    # select the particles                                             
    offset = 0
    for i in range(len(i_num)):
        for j in range(len(j_num)):
            for k in range(len(k_num)):
                box_id = divisions ** 2 * i_num[i] + divisions * j_num[j] + k_num[k]
                len_ids = len(sorted_ids[box_id])
                par_indexes[offset:offset + len_ids] = sorted_ids[box_id]
                offset += len_ids

    return par_indexes


if __name__ == "__main__":
    sims = ['m35s85', 'm4s7', 'm3s85', 'm25s9', 'm35s9', 'm35s7', 'm25s85', 'm2s8', 'm4s8', 'm2s9', 'm3s8_50', 'm3s8',
            'm35s75', 'm4s9', 'm3s9', 'm25s75', 'm2s1', 'm3s7', 'm3s75', 'm2s7', 'm25s8', 'm35s8', 'm3s8b']

    simid = int(sys.argv[1])
    mlim = float(sys.argv[2])
    snp = 118
    sim = sims[simid]  # 'm2s9'
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

    ## Load particle ids and positions

    ptype = [1]
    divisions = 200

    print('Loading pt positions and ids')
    grid_0 = np.load(folder + '1024_{}_grid_0_div{}_snap{}.npy'.format(sim, divisions, snp), allow_pickle=True)
    pos_0 = np.load(folder + '1024_{}_pos_sorted_div{}_snap{}.npy'.format(sim, divisions, snp), allow_pickle=True)
    vel_0 = np.load(folder + '1024_{}_vel_sorted_div{}_snap{}.npy'.format(sim, divisions, snp), allow_pickle=True)

    ##############-----------HALO POSITIONS----------------###########
    with open(folder + sim + '_prefixes.txt', 'r') as file:
        prefixes = file.read().splitlines()
    pdhalos = pd.read_table(folder + 'AHF/halos/{}.AHF_halos'.format(prefixes[118 - snp]), delim_whitespace=True,
                            header=0)

    print('HALO POSITIONS')

    mfof = np.array(pdhalos['Mhalo(4)'])
    hids = np.array(pdhalos['#ID(1)'])
    rvirs = np.array(pdhalos['Rhalo(12)']) / 1e3
    pos_x = np.array(pdhalos['Xc(6)']) / 1e3
    pos_y = np.array(pdhalos['Yc(7)']) / 1e3
    pos_z = np.array(pdhalos['Zc(8)']) / 1e3

    pos_h = np.array([pos_x, pos_y, pos_z]).T
    del pos_x, pos_y, pos_z

    halo_inf = []
    rmaxes = []
    nlim = 3
    for i in range(len(hids)):
        if mfof[i] > mlim:
            Rmin = 0.9 * rvirs[i]
            Rmax = 1.1 * Rmin
            Rmean = 0.5 * (Rmin + Rmax)
            Rmax_profile = nlim * Rmean
            pos_halo = pos_h[i]
            rmaxes.append(Rmax_profile)
            id_near_h0 = get_part_nearHalo(pos_halo, Rmax_profile, divisions, BoxSize, grid_0)
            pos_near_h = pos_0[id_near_h0]
            vel_near_h = vel_0[id_near_h0]
            halo_inf.append([i, pdhalos['#ID(1)'].loc[i], id_near_h0, pos_near_h, vel_near_h])

    del rvirs, hids, mfof, grid_0, pos_0, pdhalos
    if mlim == 0:
        np.save(folder + 'parts_near_hals_{}_snap{}.npy'.format(sim, snp), np.array(halo_inf, dtype=object))
    else:
        np.save(folder + 'parts_near_hals{:1.2e}_{}_snap{}.npy'.format(mlim, sim, snp),
                np.array(halo_inf, dtype=object))
