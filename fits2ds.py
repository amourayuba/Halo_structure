import numpy as np
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


#########----------------------------FUNCTIONS--------------###############
## 2D Axis ratio


def Inertia3D(points, s, q):
    '''Calculates the Inertia tensor of a set of points of equal mass in an ellipsoid where a>b>c and s=b/a and q=c/a
    points : numpy array (N,3)
    s : 0<float<=1
    q : 0<float<=1
    return 3x3 numpy array'''
    weights = points[:, 0] ** 2 + (points[:, 1] / s) ** 2 + (points[:, 2] / q) ** 2
    Iten = np.zeros((3, 3))
    for i in range(3):
        for j in range(i + 1):
            val = np.sum(points[:, i] * points[:, j] / weights)
            Iten[i, j] = val
            Iten[j, i] = val
    norm = np.trace(Iten / len(weights))  # This forces a unitary Intenria tensor
    return Iten / len(weights) / norm


def Inertia2D(points, s):
    '''Calculates the Inertian tensor of a set of points of equal mass in an ellipse where a>b and s=b/a
    points : numpy array (N,2)
    s : 0<float<=1
    return 2x2 numpy array'''
    weights = points[:, 0] ** 2 + (points[:, 1] / s) ** 2
    Iten = np.zeros((2, 2))
    for i in range(2):
        for j in range(i + 1):
            val = np.sum(points[:, i] * points[:, j] / weights)
            Iten[i, j] = val
            Iten[j, i] = val
    norm = np.trace(Iten / len(weights))
    return Iten / len(weights) / norm


def get_new_axpos2D(pos, s, selec=False):
    '''Gets the axis ratio and vector of a set of points in an ellispoid with axis a,b such as s=b/a
    pos : numpy array (N,3)
    s : 0<float <=1
    returns 
    spos : (M,3) numpy array set of positions inside the newly defined ellipse with rotated coordinates 
    news : float new value of s
    vecs : (2,2) numpy array new axis directions'''

    val, vecs = np.linalg.eig(Inertia2D(pos, s))  # Diagonalize the Inertia tensor

    arg = np.argsort(val)[::-1]  # force the first eigenvalue to be the largest for conveniance
    val = val[arg]
    vecs = vecs[:, arg]

    spos = pos.copy()

    spos[:, 0] = np.sum(pos * vecs[:, 0], axis=1)  # Rotates the x coordinates according to eigenvector 1
    spos[:, 1] = np.sum(pos * vecs[:, 1], axis=1)  # Rotates the y coordinates towards eigenvector 2

    news = np.sqrt(val[1] / val[0])

    if selec:
        spos = spos[np.sqrt(spos[:, 0] ** 2 + (spos[:, 1] / news) ** 2)]

    return spos, news, vecs


def get_axes2D(pos, selec=False, rmax=1):
    '''Gets the axes ratio s = b/a from a set of equal mass points at pos
    pos : numpy array (N,2) 
    selec : bool, whether to select only particles that are within the newly defined ellispoid. ps, not sure if selec=True works
    rmax : float, the maximim radius of particles to consider.
    Returns pos : array (N,2) final set of positions considered and rotated coordinates. 
    news : float final converged axis ratio s/q 
    vecs : numpy array (2,2) final converged axes directions'''
    # FIRST TRY
    if selec:
        pos = pos[np.sum(pos ** 2, axis=1) < rmax]
    pos, s, vecs = get_new_axpos2D(pos, 1)  # Gets the initial s startin with s = 1
    ite = 0
    while True:  # Looping over until convergence
        newpos, news, newvecs = get_new_axpos2D(pos, s, selec)  # Get the new values of s and vec
        vecs = vecs @ newvecs  # Saves the matrix by which the positions were rotated
        if (abs(news - s) / s < 1e-5):  # Condition of convergence
            return newpos, news, vecs
        else:
            # print(news)
            ite += 0
            pos, s = newpos, news
            if ite > 1000:  # Just in case it doesn't converge to get out of the loop
                print("Too many iterations, not converging.")
                return np.Nan


def plot_part_axes2D(newpart, iaxes=None, haloc=[0, 0], rvir=1):
    ''' Handy way to plot the distribution of particles and axes'''
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    ax.scatter(newpart[:, 0], newpart[:, 1], s=1)
    ax.scatter(haloc[0], haloc[1], s=10)
    c1 = plt.Circle((haloc[0], haloc[1]), rvir, color='k', fill=False)
    ax.add_patch(c1)
    ax.set_ylabel('Y [Mpc/h]')
    ax.set_xlabel('X [Mpc/h]')
    if type(iaxes) == np.ndarray:
        for j in range(2):
            xs = [-iaxes[j][0], iaxes[j][0]]
            ys = [-iaxes[j][1], iaxes[j][1]]
            ax.plot(xs, ys, color='C{}'.format(2 * j + 1), label='Axis {}'.format(j))
        ax.legend()


## $\chi^2$ of the Mass profile 

G = 4.30091e-9  # Units Mpc/Msun x (km/s)^2
rho_c = 3 * 100 ** 2 / (8 * np.pi * G)  # h^2xMsun/Mpc**3


def hubble_ratio(z, omega_l0=0.7, omega_m0=0.3, omega_r0=1e-4):
    """The value of H(z)/H0 at any given redshift for any set of present content of the
    universe"""
    omega_0 = omega_l0 + omega_m0 + omega_r0
    return np.sqrt(omega_l0 + (1 - omega_0) * (1 + z) ** 2 + omega_m0 * (1 + z) ** 3 + omega_r0 * (1 + z) ** 4)


def rhos(z, c, om0):
    '''NFW rho_s at c and z'''
    deltac = 200 * c ** 3 / (3 * (np.log(1 + c) - c / (1 + c)))
    return deltac * rho_c * hubble_ratio(z, omega_l0=1 - om0, omega_m0=om0) ** 2


def rhotilde(r, loga):
    return 1 / ((np.log(2) - 0.5) * (r / 10 ** loga) * (1 + r / 10 ** loga) ** 2)


def Ntilde(r, loga):
    return (np.log(1 + r / 10 ** loga) - (r / 10 ** loga) / (1 + r / 10 ** loga)) / (np.log(2) - 0.5)


def prob3D(r, rmin, rmax, loga):
    return (r / 10 ** loga) ** 2 * rhotilde(r, loga) / (10 ** loga * (Ntilde(rmax, loga) - Ntilde(rmin, loga)))


def sigma_tilde(R, loga):
    """Normalized surface density function for an NFW profile. 
    loga: float, scale radius of the wanted 2D NFW profile. 
    R: positive float or array of positive floats, Distance from the center of the profile
    Returns: float, or array of floats the value of the normalized surface density profile"""

    nor = 0.5 / (1 - np.log(2))  # using a temporary variable to make the math clear
    a, X = 10 ** loga, R / 10 ** loga

    if type(X) == np.ndarray or type(X) == list:  # case R is an array, to get an array as result
        X1 = X[np.where(X > 1)]
        resu3 = nor * (1 - np.arccos(1 / X1) / np.sqrt(X1 ** 2 - 1)) / (X1 ** 2 - 1)
        X2 = X[np.where((X < 1) & (X > 0))]
        resu1 = nor * (1 - np.arccosh(1 / X2) / np.sqrt(-X2 ** 2 + 1)) / (X2 ** 2 - 1)

        if len(X[np.where(X == 1)]) == 1:  # Case Ri = a
            resu2 = np.array([nor * 1 / 3])
            return np.concatenate((resu1, resu2, resu3))
        return np.concatenate((resu1, resu3))
    # if R is a unique number

    if X > 1:
        return nor * (1 - np.arccos(1 / X) / np.sqrt(X ** 2 - 1)) / (X ** 2 - 1)
    elif X < 1:
        return nor * (1 - np.arccosh(1 / X) / np.sqrt(-X ** 2 + 1)) / (X ** 2 - 1)
    elif X == 1:
        return nor * 1 / 3
    else:
        print("R has to be positive")


def N_tilde2D(R, loga):
    """Normalized Number of particles function for an NFW profile. 
    loga: float, scale radius of the wanted 2D NFW profile. 
    R: positive float or array of positive floats, Distance from the center of the profile
    Returns: float, or array of floats the value of the normalized surface density profile"""

    nor = 1 / (1 - np.log(2))  # using a temporary variable to make the math clear
    a, X = 10 ** loga, R / 10 ** loga

    if type(X) == np.ndarray or type(X) == list:

        l0 = len(np.where(X == 0)[0])
        resu0 = [0] * l0

        X1 = X[np.where(X > 1)]
        resu3 = nor * (np.arccos(1 / X1) / np.sqrt(X1 ** 2 - 1) + np.log(X1 / 2))

        X2 = X[np.where((X < 1) & (X > 0))]
        resu1 = nor * (np.arccosh(1 / X2) / np.sqrt(-X2 ** 2 + 1) + np.log(X2 / 2))

        if len(X[np.where(X == 1)]) == 1:
            resu2 = np.array([nor * (1 - np.log(2))])
            return np.concatenate((resu1, resu2, resu3))
        return np.concatenate((resu0, resu1, resu3))
    if X == 0:
        return 0
    if X > 1:
        return nor * (np.arccos(1 / X) / np.sqrt(X ** 2 - 1) + np.log(X / 2))
    elif X < 1:
        return nor * (np.arccosh(1 / X) / np.sqrt(-X ** 2 + 1) + np.log(X / 2))
    elif X == 1:
        return nor * (1 - np.log(2))
    else:
        print("R has to be positive")


def NFW3d(R, c, rvir, z, om0):
    norm = rhos(z, c, om0)
    rs = rvir / c
    return norm / ((R / rs) * (1 + R / rs) ** 2)


def Mass_NFW(r, c, Mvir, rvir):
    rs = rvir / c
    return Mvir * (np.log(1 + r / rs) - (r / rs) / (1 + r / rs)) / (np.log(1 + c) - c / (1 + c))


def Mass_NFW2d(R, c, rvir, om0, z=0):
    rs = rvir / c
    norm = 4 * np.pi * rs ** 3 * rhos(z, c, om0) * (1 - np.log(2))
    return norm * N_tilde2D(R, np.log10(rs))


def NFW2d(R, c, rvir, om0, z=0):
    rs = rvir / c
    norm = 2 * rs * rhos(z, c, om0) * (2 - 2 * np.log(2))
    X = R / rs
    return norm * sigma_tilde(R, np.log10(rs))


def chi_square(f, x_data, y_data, unc):
    '''Gives the chi**2 of a set of points following the function f with the wiki def of chi**2'''
    # return np.sum((f(x_data)-y_data)**2/f(x_data))
    return np.sum((f(x_data) - y_data) ** 2 / unc ** 2)


def log_chi_square(f, x_data, y_data):
    '''Gives the chi**2 of a set of points following the function f with the more reasonable def of chi**2'''
    nbins = len(x_data)
    return np.sum((np.log10(f(x_data) / y_data)) ** 2) / nbins


## Concentration minimisation

def prob2D(R, Rmin, Rmax, loga):
    # print(loga, len(R), loga.shape)
    a1 = 2 * (R / 10 ** loga)
    # print(a1.shape)
    a2 = sigma_tilde(R, loga)
    # print(a2.shape)
    a3 = N_tilde2D(Rmax, loga) - N_tilde2D(Rmin, loga)
    # print(loga, Rmax, N_tilde2D(Rmax, loga), Rmin, N_tilde2D(Rmin, loga))
    a3b = 10 ** loga * a3
    # print(a3b.shape)
    return a1 * a2 / a3b


def N_frac(R, dmax, a):
    return N_tilde(R, np.log10(a)) / N_tilde(dmax, np.log10(a))


def likelihood(rmin, rmax, loga, func, rad_arrays):
    likelihood = 0
    indices = np.where((rad_arrays >= rmin) & (rad_arrays <= rmax))
    likelihood = np.sum(-np.log(func(rad_arrays[indices], rmin, rmax, loga)))
    # print (len(indices[0]))
    return likelihood


def minimisation(fonction, rmin, rmax, func, rad_arrays, x0, argbounds=None, methode=None):
    def f(x):
        return fonction(rmin, rmax, x, func, rad_arrays)

    if methode == None:
        return minimize(f, x0, bounds=argbounds).x
    else:
        return minimize(f, x0, method=methode, bounds=argbounds).x

    ## Most bound particle


def g_pot_r(particle_pos, halo_c, nbins=1000, mpart=9690706183.9515):
    '''Calculates the gravitation potential of equal mass set of particles'''

    prads = np.sqrt(np.sum((particle_pos - halo_c) ** 2, axis=1))  # calculates the radius

    rbounds = np.logspace(-3, -1, nbins)  # min is softening length, max is arbitrary, but too far is irrelevant

    h, rb = np.histogram(prads, bins=rbounds)
    Mass_inr = np.cumsum(h) * mpart  # M(<r)
    pot = G * np.cumsum(Mass_inr * (rb[1:] - rb[:-1]) / rb[1:] ** 2)  # The G potential
    return prads, rb, pot


def get_binding_e(p_pos, vel_sel, halo_c, nbins=1000, mpart=9690706183.9515):
    '''Calculates the binding Energy (not considering phi0, total halo Energy)'''
    rad_sel, rb, pot = g_pot_r(p_pos, halo_c, nbins, mpart)  # Potential Phi(r)
    argbin = np.digitize(rad_sel, rb[1:-1])  # Gives where each particle's radial bin is

    part_pot = pot[argbin]  # Gives each particle's potential

    return rad_sel, np.sum(vel_sel ** 2, axis=1) / 2 + part_pot


##################---------------GETTING THE PROPERTY-------############################## 


def load_parts(halo_inf, pdhalos, halid, rvirlim=3):
    halo = halo_inf[halid]
    ahfid = halo[0]
    vels, parts = halo[-1].astype(np.float64), halo[-2].astype(np.float64)

    xh, yh, zh, rvir = pdhalos['Xc(6)'].loc[ahfid] / 1e3, pdhalos['Yc(7)'].loc[ahfid] / 1e3, pdhalos['Zc(8)'].loc[
        ahfid] / 1e3, pdhalos['Rhalo(12)'].loc[ahfid] / 1e3
    vxh, vyh, vzh, mvir = pdhalos['VXc(9)'].loc[ahfid], pdhalos['VYc(10)'].loc[ahfid], pdhalos['VZc(11)'].loc[ahfid], \
    pdhalos['Mhalo(4)'].loc[ahfid]

    part_pos = parts.copy()
    if xh % 500 < rvirlim * rvir:
        xh = (xh - 250) % 500
        part_pos[:, 0] = (parts[:, 0] - 250) % 500
    if yh % 500 < rvirlim * rvir:
        yh = (yh - 250) % 500
        part_pos[:, 1] = (parts[:, 1] - 250) % 500
    if zh % 500 < rvirlim * rvir:
        zh = (zh - 250) % 500
        part_pos[:, 2] = (parts[:, 2] - 250) % 500

    rads3d = np.sqrt(np.sum((part_pos - np.array([xh, yh, zh])) ** 2, axis=1))
    rads2d = [np.sqrt((part_pos[:, 0] - xh) ** 2 + (part_pos[:, 1] - yh) ** 2),
              np.sqrt((part_pos[:, 0] - xh) ** 2 + (part_pos[:, 2] - zh) ** 2),
              np.sqrt((part_pos[:, 1] - yh) ** 2 + (part_pos[:, 2] - zh) ** 2)]

    return vels, part_pos, np.array([xh, yh, zh]), np.array([vxh, vyh, vzh]), rvir, mvir, rads3d, rads2d


def get_axe(part_pos, halo_c, rvir):
    newpart = (part_pos[:, :-1] - halo_c[:-1]) / rvir
    newpart2 = (part_pos[:, 1:] - halo_c[1:]) / rvir
    newpart3 = (part_pos[:, ::2] - halo_c[::2]) / rvir

    parts_xy = newpart[np.sum(newpart ** 2, axis=1) <= 1]
    parts_yz = newpart[np.sum(newpart2 ** 2, axis=1) <= 1]
    parts_xz = newpart[np.sum(newpart3 ** 2, axis=1) <= 1]

    newpos, news, vecs = get_axes2D(parts_xy)
    newpos2, news2, vecs2 = get_axes2D(parts_yz)
    newpos3, news3, vecs3 = get_axes2D(parts_xz)

    return np.array([news, news2, news3]), np.array([vecs, vecs2, vecs3])


def get_mbp_off1(part_pos, part_vel, halo_c, vel_h):
    pr, bE = get_binding_e(part_pos, part_vel - vel_h, halo_c, nbins=1000)
    xy_R = np.sqrt(np.sum((part_pos[:, :2] - halo_c[:2]) ** 2, axis=1))
    yz_R = np.sqrt(np.sum((part_pos[:, 1:] - halo_c[1:]) ** 2, axis=1))
    xz_R = np.sqrt(np.sum((part_pos[:, ::2] - halo_c[::2]) ** 2, axis=1))

    mbp = np.argmin(bE)
    return np.array([xy_R[mbp], yz_R[mbp], xz_R[mbp]]), pr[mbp]


def get_com_off1(part_pos, rads, rvir, halo_c):
    selec_pos = np.where(rads <= rvir)
    spos = part_pos[selec_pos]
    com_pos = np.array([np.mean(spos[:, 0]), np.mean(spos[:, 1]), np.mean(spos[:, 2])])
    com3d = np.sqrt(np.sum((com_pos - halo_c) ** 2))

    comxy = np.sqrt(np.sum((com_pos[:2] - halo_c[:2]) ** 2))
    comyz = np.sqrt(np.sum((com_pos[1:] - halo_c[1:]) ** 2))
    comxz = np.sqrt(np.sum((com_pos[::2] - halo_c[::2]) ** 2))
    return com3d, np.array([comxy, comyz, comxz])


def get_conc3d(radii, rvir, rmin=0, rmax=1.2):
    resu = \
    minimisation(likelihood, rmin * rvir, rmax * rvir, prob3D, radii, np.log10(np.median(radii)), methode='COBYLA')[0]
    return rvir / 10 ** resu


def get_conc2d(rads, rvir, rmin=0, rmax=1.2):
    res = []
    for radii in rads:
        resu = \
        minimisation(likelihood, rmin * rvir, rmax * rvir, prob2D, radii, np.log10(np.median(radii)), methode='COBYLA')[
            0]
        res.append(rvir / 10 ** resu)
    return np.array(res)


def get_chi3d(conc, rads, rvir, mvir, z=0, nbins=10, mpart=9690706183.9515, om0=0.3):
    redbins = np.logspace(-2, np.log10(2 * rvir), nbins)
    binvolume = 4 * np.pi * redbins[1:] * redbins[:-1] * (redbins[1:] - redbins[:-1])  # Mpc**3
    hist, rbins = np.histogram(rads, bins=redbins)
    rho_d = mpart * hist / binvolume
    rho_u = mpart * np.sqrt(hist) / binvolume
    mass_p = mpart * np.cumsum(hist)
    mass_u = mpart * np.sqrt(np.cumsum(hist))
    c_redbin = np.sqrt(redbins[1:] * redbins[:-1])

    def fm(x):
        return Mass_NFW(x, conc, mvir, rvir)

    def frho(x):
        return NFW3d(x, conc, rvir, z, om0)

    chim = np.array([chi_square(fm, c_redbin, mass_p, mass_u) / mvir, log_chi_square(fm, c_redbin, mass_p)])
    chirho = np.array([chi_square(frho, c_redbin, rho_d, rho_u) / mvir, log_chi_square(frho, c_redbin, rho_d)])
    return chim, chirho


def get_chi2d(conc, rads, rvir, mvir, z=0, nbins=10, mpart=9690706183.9515, om0=0.3):
    redbins = np.logspace(-2, np.log10(2 * rvir), nbins)
    binsurf = 2 * np.pi * np.sqrt(redbins[1:] * redbins[:-1]) * (redbins[1:] - redbins[:-1])  # Mpc**3
    hist, rbins = np.histogram(rads, bins=redbins)

    sig_d = mpart * hist / binsurf
    sig_u = mpart * np.sqrt(hist) / binsurf
    mass_p = mpart * np.cumsum(hist)
    mass_u = mpart * np.sqrt(np.cumsum(hist))
    c_redbin = np.sqrt(redbins[1:] * redbins[:-1])

    def fm(x):
        return Mass_NFW2d(x, conc, rvir, z, om0)

    def frho(x):
        return NFW2d(x, conc, rvir, z, om0)

    chim = np.array([chi_square(fm, c_redbin, mass_p, mass_u) / mvir, log_chi_square(fm, c_redbin, sig_d)])

    chirho = np.array([chi_square(frho, c_redbin, sig_d, sig_u) / mvir, log_chi_square(frho, c_redbin, sig_d)])
    return chim, chirho
