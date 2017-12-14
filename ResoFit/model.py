import numpy as np
import matplotlib.pyplot as plt
from lmfit.lineshapes import pvoigt


def ikeda_carpenter(t, alpha, beta, fraction, t0):
    """

    :param t:
    :type t:
    :param alpha:
    :type alpha:
    :param beta:
    :type beta:
    :param t0:
    :type t0:
    :param fraction:
    :type fraction:
    :return:
    :rtype:
    """
    _t = t - t0
    part1 = (1 - fraction) * (alpha * _t) ** 2 * np.exp(-alpha * _t)
    part2 = np.exp(-beta * _t) - np.exp(-alpha * _t) * (
            1 + (alpha - beta) * _t + 0.5 * ((alpha - beta) ** 2) * (_t ** 2))

    f = 0.5 * alpha * (part1 + (2 * fraction * alpha * beta / (alpha - beta) ** 3) * part2)
    return f


def cole_windsor(t, sig1, sig2, gam1, gam2, norm_factor, fraction, t0):
    """


    :param t:
    :type t:
    :param sig1:
    :type sig1:
    :param sig2:
    :type sig2:
    :param gam1:
    :type gam1:
    :param gam2:
    :type gam2:
    :param norm_factor:
    :type norm_factor:
    :param fraction:
    :type fraction:
    :param t0:
    :type t0:
    :return:
    :rtype:
    """
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

# def pseudo_voigt()


x = np.linspace(0, 20, 1000)
# y = ikeda_carpenter(x, 0.5, 1.3, 1.7, 1.05)
y = cole_windsor(x, 0.5, 1.3, 1.7, 1.05, 1.7, 1, 2)

plt.plot(x, y)
plt.show()
