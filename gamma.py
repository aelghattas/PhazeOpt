import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma
from scipy.optimize import curve_fit
plt.style.use('tableau-colorblind10')

class Gamma:
    def __init__(self, filename, sheet, scn_cell, mw_cell, mole_fraction_cell, h_parameter=-4):
        self.h = 0
        self.scn_start = 0
        self.scn_plus_start = 0
        self.alpha = 0
        self.beta = 0
        self.lower_bound = 0
        self.mi_char_lab = 0
        self.zni_char_lab = 0
        self.zn_plus_lab = 0
        self.mw_plus_lab = 0
        self.zn_plus_split = 0
        self.mw_plus_split = 0
        self.maximum_split = 500
        self.split_fit_mi_full = 0
        self.split_fit_zni_full = 0
        self.data = pd.read_excel(filename, sheet_name=sheet)
        self.scn = np.array(self.data[scn_cell])
        self.mi_lab = np.array(self.data[mw_cell])
        self.mf_lab = np.array(self.data[mole_fraction_cell])
        self.h = h_parameter
        self.characterized(self.scn, self.mi_lab, self.mf_lab, self.h)


    def gamma_parameters(self, x, alpha_, beta_):
        return gamma.cdf(x / beta_, alpha_)
    
    def split(self, scn_start_, scn_plus_start_, scn_plus_end_, h_, alpha_, beta_, lower_bound_):
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

    def split_error(self, split_limit):
        output = self.split(self.scn_start, self.scn_plus_start, split_limit, self.h, self.alpha, self.beta, self.lower_bound)
        return (output[0] - self.zn_plus_lab)**2


    def optimize_split(self, maximum_split_):
        optimum_ = 0
        max_split_range = np.arange(self.scn_plus_start + 3, maximum_split_)
        for i, j in enumerate(max_split_range):
            if i > 0 and self.split_error(j) - self.split_error(j - 1) > 0:
                optimum_ = (j - 1)
                break
            else:
                continue
        return optimum_
    
    def characterized(self, scn_in, mi_in, mf_in, h_in):
        zni_lab = mf_in / mf_in.sum()
        zni_cum_lab = np.cumsum(zni_lab)

        self.zni_char_lab = zni_lab[:-2]
        self.mi_char_lab = mi_in[:-2]

        self.zn_plus_lab = zni_lab[-2:].sum()
        self.mw_plus_lab = np.dot(mi_in[-2:], zni_lab[-2:]) / self.zn_plus_lab

        self.scn_start = scn_in[0]
        self.scn_plus_start = scn_in[-2]

        self.lower_bound = np.array(scn_in[0] * 14) + h_in - 7
        x_cdf_char = np.array(scn_in * 14) + h_in + 7 - self.lower_bound

        fit_parameters = curve_fit(self.gamma_parameters, x_cdf_char, zni_cum_lab)
        self.alpha = fit_parameters[0][0]
        self.beta = fit_parameters[0][1]

        self.maximum_split = self.optimize_split(self.maximum_split)
        self.zn_plus_split, self.mw_plus_split = self.split(self.scn_start, self.scn_plus_start, self.maximum_split, self.h, self.alpha, self.beta, self.lower_bound)

        scn_split_full = np.arange(self.scn_start, self.maximum_split + 1)
        self.split_fit_mi_full = np.array(scn_split_full * 14) + h_in
        x_cdf_fit = self.split_fit_mi_full + 7 - self.lower_bound
        zni_fit_cdf = gamma.cdf(x_cdf_fit / self.beta, self.alpha)
        zni_fit_cdf[1:] -= zni_fit_cdf[:-1].copy()
        self.split_fit_zni_full = zni_fit_cdf

    def report(self):
        print('Optimized Parameters: ')
        print('>> Optimized Split SCN: ' + str(self.maximum_split))
        print('>> Bound: ' + str(self.lower_bound))
        print('>> Alpha: ' + str(self.alpha))
        print('>> Beta: ' + str(self.beta))
        print('>> h value used for Mi Calculations: ' + str(self.h))

    def plot(self):
        plt.scatter(self.mi_char_lab, self.zni_char_lab, label='Lab Data')
        plt.plot(self.split_fit_mi_full, self.split_fit_zni_full, label='Fit Data')
        plt.scatter(self.mw_plus_lab, self.zn_plus_lab, label='Lab Plus Fraction')
        plt.scatter(self.mw_plus_split, self.zn_plus_split, label='Fit Plus Fraction')
        plt.xlabel('Molecular Weight')
        plt.ylabel('Normalized Mole Fraction')
        plt.xlim(0)
        plt.ylim(0)
        plt.legend()
        plt.show()        
