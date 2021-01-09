import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma
from scipy.optimize import curve_fit
plt.style.use('tableau-colorblind10')

# global variables
h = 0
scn_start = 0
scn_plus_start = 0
alpha = 0
beta = 0
lower_bound = 0
mi_char_lab = 0
zni_char_lab = 0
zn_plus_lab = 0
mw_plus_lab = 0
zn_plus_split = 0
mw_plus_split = 0
maximum_split = 500
split_fit_mi_full = 0
split_fit_zni_full = 0


def gamma_parameters(x, alpha_, beta_):
    return gamma.cdf(x / beta_, alpha_)


def split(scn_start_, scn_plus_start_, scn_plus_end_, h_, alpha_, beta_, lower_bound_):
    scn_split = np.arange(scn_start_, scn_plus_end_)
    mi_split = np.array(scn_split * 14) + h_
    x_cdf_split = mi_split + 7 - lower_bound_
    zni_split_cdf = gamma.cdf(x_cdf_split / beta_, alpha_)
    _zn_split = zni_split_cdf
    _zn_split[1:] -= _zn_split[:-1].copy()
    _zn_plus_split = _zn_split[(scn_plus_start_ - scn_start_):].sum()
    _mw_plus_split = np.dot(mi_split[(scn_plus_start_ - scn_start_):], _zn_split[(scn_plus_start_ - scn_start_):]) / \
        _zn_plus_split
    return _zn_plus_split, _mw_plus_split


def split_error(split_limit):
    output = split(scn_start, scn_plus_start, split_limit, h, alpha, beta, lower_bound)
    return (output[0] - zn_plus_lab)**2


def optimize_split(maximum_split_):
    optimum_ = 0
    max_split_range = np.arange(scn_plus_start + 3, maximum_split_)
    for i, j in enumerate(max_split_range):
        if i > 0 and split_error(j) - split_error(j - 1) > 0:
            optimum_ = (j - 1)
            break
        else:
            continue
    return optimum_


def characterized(scn_in, mi_in, mf_in, h_in):
    global scn_start
    global scn_plus_start
    global alpha
    global beta
    global lower_bound
    global mi_char_lab
    global zni_char_lab
    global zn_plus_lab
    global mw_plus_lab
    global zn_plus_split
    global mw_plus_split
    global maximum_split
    global split_fit_mi_full
    global split_fit_zni_full
    global h

    h = h_in

    zni_lab = mf_in / mf_in.sum()
    zni_cum_lab = np.cumsum(zni_lab)

    zni_char_lab = zni_lab[:-2]
    mi_char_lab = mi_in[:-2]

    zn_plus_lab = zni_lab[-2:].sum()
    mw_plus_lab = np.dot(mi_in[-2:], zni_lab[-2:]) / zn_plus_lab

    scn_start = scn_in[0]
    scn_plus_start = scn_in[-2]

    lower_bound = np.array(scn_in[0] * 14) + h_in - 7
    x_cdf_char = np.array(scn_in * 14) + h_in + 7 - lower_bound

    fit_parameters = curve_fit(gamma_parameters, x_cdf_char, zni_cum_lab)
    alpha = fit_parameters[0][0]
    beta = fit_parameters[0][1]


    maximum_split = optimize_split(maximum_split)
    zn_plus_split, mw_plus_split = split(scn_start, scn_plus_start, maximum_split, h, alpha, beta, lower_bound)

    scn_split_full = np.arange(scn_start, maximum_split + 1)
    split_fit_mi_full = np.array(scn_split_full * 14) + h_in
    x_cdf_fit = split_fit_mi_full + 7 - lower_bound
    zni_fit_cdf = gamma.cdf(x_cdf_fit / beta, alpha)
    zni_fit_cdf[1:] -= zni_fit_cdf[:-1].copy()
    split_fit_zni_full = zni_fit_cdf


def report():
    print('Optimized Parameters: ')
    print('>> Optimized Split SCN: ' + str(maximum_split))
    print('>> Bound: ' + str(lower_bound))
    print('>> Alpha: ' + str(alpha))
    print('>> Beta: ' + str(beta))
    print('>> h value used for Mi Calculations: ' + str(h))


def fit(filename, sheet, scn_cell, mw_cell, mole_fraction_cell, h_parameter=-4):
    data = pd.read_excel(filename, sheet_name=sheet)
    scn = np.array(data[scn_cell])
    mi_lab = np.array(data[mw_cell])
    mf_lab = np.array(data[mole_fraction_cell])
    characterized(scn, mi_lab, mf_lab, h_parameter)
    report()
    plt.scatter(mi_char_lab, zni_char_lab, label='Lab Data')
    plt.plot(split_fit_mi_full, split_fit_zni_full, label='Fit Data')
    plt.scatter(mw_plus_lab, zn_plus_lab, label='Lab Plus Fraction')
    plt.scatter(mw_plus_split, zn_plus_split, label='Fit Plus Fraction')
    plt.xlabel('Molecular Weight')
    plt.ylabel('Normalized Mole Fraction')
    plt.xlim(0)
    plt.ylim(0)
    plt.legend()
    plt.show()


if __name__ == '__main__':
    fit('data.xlsx', 'Data', 'SCN', 'Lab Mi', 'mol-%', -6)
