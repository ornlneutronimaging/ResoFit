import numpy as np
import matplotlib.pyplot as plt
from lmfit.lineshapes import pvoigt
from lmfit import Model


def ikeda_carpenter(t, alpha, beta, fraction, t0, norm_factor=1):
    _t = t - t0
    # _t = t1[np.logical_not(t1 < 0)]
    _t[_t < 0] = 0  # t>=0 required

    # α=vΣ
    # Σ is the macroscopic neutron scattering cross-section of the moderator
    #   (Σ=1.6 cm-1 for polyethylene in the high energy limit) and
    # v is the neutron speed

    part1 = 0.5 * alpha * (alpha*_t)**2 * np.exp(-alpha*_t)
    part2_1 = alpha**3 * beta / (alpha - beta)**3  # jparc uses alpha**2 instead of alpha**3 in the original function
    part2_2 = np.exp(-beta*_t) - np.exp(-alpha*_t) * (1 + (alpha - beta) * _t + 0.5 * (alpha - beta)**2 * _t**2)
    part2 = part2_1 * part2_2

    f = ((1 - fraction) * part1 + fraction * part2) * norm_factor
    return f


def ikeda_carpenter_jparc(t, alpha, beta, fraction, t0, norm_factor=1):
    _t = t - t0
    # _t = t1[np.logical_not(t1 < 0)]
    _t[_t < 0] = 0  # t>=0 required

    part1 = 0.5 * alpha * (alpha*_t)**2 * np.exp(-alpha*_t)
    part2_1 = alpha**2 * beta / (alpha - beta)**3  # jparc uses alpha**2 instead of alpha**3 in the original function
    part2_2 = np.exp(-beta*_t) - np.exp(-alpha*_t) * (1 + (alpha - beta) * _t + 0.5 * (alpha - beta)**2 * _t**2)
    part2 = part2_1 * part2_2

    f = ((1 - fraction) * part1 + fraction * part2) * norm_factor
    return f


def cole_windsor(t, sig1, sig2, gamma, fraction, t0, norm_factor=1):
    _t = t - t0
    f = []
    for each_t in _t:
        # for F1
        if each_t <= 0:
            f1 = np.exp(-0.5 * (each_t / sig1) ** 2)
        if 0 < each_t <= gamma * sig2 ** 2:
            f1 = np.exp(-0.5 * (each_t / sig2) ** 2)
        if gamma * sig2 ** 2 < each_t:
            f1 = np.exp(0.5 * (gamma * sig2) ** 2 - gamma * each_t)
        # for F2
        if each_t <= 0:
            f2 = np.exp(-0.5 * (each_t / sig1) ** 2)
        if 0 < each_t <= gamma * sig2 ** 2:
            f2 = np.exp(0.5 * (each_t / sig2) ** 2)
        if gamma * sig2 ** 2 < each_t:
            f2 = np.exp(0.5 * (gamma * sig2) ** 2 - gamma * each_t)

        each_f = norm_factor * ((1 - fraction) * f1 + fraction * f2)
        f.append(each_f)

    f = np.array(f)
    return f


def pseudo_voigt(t, beta, sigma, fraction):
    gauss = 1 / (1 + (t / beta)**2)
    lorentz = np.exp(-(t / sigma)**2)
    f = (1 - fraction) * gauss + fraction * lorentz
    return f


def cole_windsor_jparc(t, sig1, sig2, gam1, gam2, fraction, t0, norm_factor=1):
    _t = t - t0
    f = []
    for each_t in _t:
        # for F1
        if each_t <= 0:
            f1 = np.exp(-0.5 * (each_t / sig1) ** 2)
        if 0 < each_t <= gam1 * sig2 ** 2:
            f1 = np.exp(-0.5 * (each_t / sig2) ** 2)
        if gam1 * sig2 ** 2 < each_t:
            f1 = np.exp(0.5 * (gam1 * sig2) ** 2 - gam1 * each_t)
        # for F2
        if each_t <= 0:
            f2 = np.exp(-0.5 * (each_t / sig1) ** 2)
        if 0 < each_t <= gam2 * sig2 ** 2:
            f2 = np.exp(0.5 * (each_t / sig2) ** 2)
        if gam2 * sig2 ** 2 < each_t:
            f2 = np.exp(0.5 * (gam2 * sig2) ** 2 - gam2 * each_t)

        each_f = norm_factor * ((1 - fraction) * f1 + fraction * f2)
        f.append(each_f)

    f = np.array(f)
    return f


