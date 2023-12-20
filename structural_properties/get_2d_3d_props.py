import numpy as np
import pandas as pd
import sys
import fits2ds as f2d


if __name__ == "__main__":
    sims = ['m4s7', 'm3s85', 'm25s9', 'm35s9', 'm35s7', 'm25s85', 'm2s8', 'm4s8', 'm2s9', 'm3s8', 'm35s75', 'm4s9', 'm3s9', 'm25s75', 'm2s1', 'm3s7', 'm3s75', 'm2s7', 'm25s8', 'm35s8', 'm35s85']
    oms = [0.4, 0.3, 0.25, 0.35, 0.35, 0.25, 0.2, 0.4, 0.2, 0.3, 0.35, 0.4, 0.3, 0.25, 0.2, 0.3, 0.3, 0.2, 0.25, 0.35, 0.35]

    simid = int(sys.argv[1])
    mlim = float(sys.argv[2])
    snp = 118
    sim = sims[simid] #'m2s9'
    folder = '/home/ayuba/scratch/{}/'.format(sim)
    print(folder)
    ## Load particle ids and positions


    with open(folder+sim+'_prefixes.txt', 'r') as file:
        prefixes = file.read().splitlines()
    pdhalos = pd.read_table(folder+'AHF/halos/{}.AHF_halos'.format(prefixes[118-snp]), delim_whitespace=True, header=0)

    mfof = np.array(pdhalos['Mhalo(4)'])
    hids = np.array(pdhalos['#ID(1)'])

    if mlim == 0:
        halo_inf = np.load(folder+'parts_near_hals_{}_snap{}.npy'.format(sim, snp), allow_pickle=True)
    else:
        halo_inf = np.load(folder+'parts_near_hals{:1.2e}_{}_snap{}.npy'.format(mlim, sim, snp))


    props3d = []
    props2d = []
    k = 0
    for i in range(len(hids)):
        if mfof[i] > mlim:
            vels, part_pos, halo_c, vel_h, rvir, mvir, rads3d, rads2d = f2d.load_parts(halo_inf, pdhalos, i)
            axs2d, vecs2d = f2d.get_axe(part_pos, halo_c, rvir)
            mboff2d, _ = f2d.get_mbp_off1(part_pos, vels, halo_c, vel_h)
            com3d, com2d = f2d.get_com_off1(part_pos, rads3d, rvir, halo_c)
            conc3d = f2d.get_conc3d(rads3d, rvir, rmin=0, rmax=1.2)
            conc2d = f2d.get_conc2d(rads2d, rvir, rmin=0, rmax=1.2)
            chim3d, chirho3d = f2d.get_chi3d(conc3d, rads3d, rvir, mvir, z=0, nbins=10, om0=oms[simid])
            chim2d, chirho2d = f2d.get_chi2d(conc2d[k], rads2d[k], rvir, mvir, z=0, nbins=10, om0=oms[simid]) #xyplane

            props2d.append([conc2d[k], chirho2d, chim2d, axs2d[k], vecs2d[k], mboff2d[k], com2d[k]])
            props3d.append([conc3d, chirho3d, chim3d])

            if i%1000 == 0:
                print('Halo number ', i)
    props3d = np.array(props3d, dtype=object)
    props2d = np.array(props2d, dtype=object)

    if mlim == 0:
        np.save(folder+'nprops3D_{}_snap{}.npy'.format(sim, snp), props3d)
        np.save(folder+'nprops2D_{}_snap{}.npy'.format(sim, snp), props2d)
    else:
        np.save(folder+'nprops3D{:1.2e}_{}_snap{}.npy'.format(mlim, sim, snp), props3d)
        np.save(folder+'nprops2D{:1.2e}_{}_snap{}.npy'.format(mlim, sim, snp), props2d)
