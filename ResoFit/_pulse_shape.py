import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import Parameters
from lmfit import minimize
from ResoFit._gap_functions import gap_neutron_pulse_ikeda_carpenter
from ResoFit._gap_functions import gap_neutron_pulse_cole_windsor
import ImagingReso._utilities as reso_util
import ResoFit._utilities as fit_util


class NeutronPulse(object):

    def __init__(self, path):
        """"""
        self.shape_total_df = load_neutron_total_shape(path)
        self.params_to_fitshape = None
        self.shape_result = None
        self.shape_dict = None
        self.model_index = None

    def load_shape_each(self, path):
        self.shape_dict = load_neutron_each_shape(path)

    def export_total(self, filename=None):
        assert self.shape_total_df is not None

        if filename is None:
            self.shape_total_df.to_clipboard(excel=True)
        else:
            self.shape_total_df.to_csv(filename)

    def export_each(self):
        for index, each_energy in enumerate(self.shape_dict.keys()):
            df = self.shape_dict[each_energy]
            file_name = 'energy_' + str(index + 1) + '.csv'
            df.to_csv(file_name, index=False)
            print("Neutron pulse shape of 'E = {} eV' has exported to './{}'".format(each_energy, file_name))

    def fit_shape(self, t, f, model_index=1, each_step=False):
        # [1: 'ikeda_carpenter', 2: 'cole_windsor', 3: 'pseudo_voigt']
        self.model_index = model_index
        self.params_to_fitshape = Parameters()

        # ikeda_carpenter
        if self.model_index == 1:
            # Load params
            self.params_to_fitshape.add('alpha', value=0.06)
            self.params_to_fitshape.add('beta', value=0.05)
            self.params_to_fitshape.add('fraction', value=0.5, min=0, max=1)
            self.params_to_fitshape.add('t0', value=0.01)
            # Use lmfit to obtain params by minimizing gap_function
            self.shape_result = minimize(gap_neutron_pulse_ikeda_carpenter,
                                         self.params_to_fitshape,
                                         method='leastsq',
                                         args=(t, f, each_step)
                                         )
        # cole_windsor
        elif self.model_index == 2:
            # Load params
            self.params_to_fitshape.add('sig1',
                                        # value=source_to_detector_m
                                        )
            self.params_to_fitshape.add('sig2',
                                        # value=offset_us
                                        )
            self.params_to_fitshape.add('gam1',
                                        # value=offset_us
                                        )
            self.params_to_fitshape.add('gam2',
                                        # value=offset_us
                                        )
            self.params_to_fitshape.add('norm_factor',
                                        # value=source_to_detector_m
                                        )
            self.params_to_fitshape.add('fraction',
                                        # value=0.5,
                                        min=0,
                                        max=1
                                        )
            self.params_to_fitshape.add('t0',
                                        # value=offset_us,
                                        vary=True
                                        )
            # Use lmfit to obtain params by minimizing gap_function
            self.shape_result = minimize(gap_neutron_pulse_cole_windsor,
                                         self.params_to_fitshape,
                                         method='leastsq',
                                         args=(t, f, each_step))

        # Print before
        print("+----------------- Fit neutron pulse shape -----------------+\nParams before:")
        self.params_to_fitshape.pretty_print()
        # Use lmfit to obtain params by minimizing gap_function

        # Print after
        print("\nParams after:")
        self.shape_result.__dict__['params'].pretty_print()
        # Print chi^2
        print("Calibration chi^2 : {}\n".format(self.shape_result.__dict__['chisqr']))


class ProtonPulse(object):
    pass


def load_neutron_total_shape(path):
    p = np.genfromtxt(path)
    pp = p.T
    df1 = pd.DataFrame()
    for i, col in enumerate(pp):
        df1[i] = col
    col_name_1 = ['E_eV', 'I_angstrom', 'f(E)', 's(E)', 'f(I)', 's(I)']
    df1.columns = col_name_1
    return df1


def load_neutron_each_shape(path):
    q = np.genfromtxt(path)
    qq = q.T
    df2 = pd.DataFrame()
    for i, col in enumerate(qq):
        df2[i] = col

    energy_list = list(set(list(df2[1])))
    energy_list.sort()

    shape_dict = {}
    for index, each_energy in enumerate(energy_list):
        t_us = []
        e_ev = []
        f = []
        s = []
        df = pd.DataFrame()
        for each_line in q:
            if each_energy == each_line[1]:
                t_us.append(each_line[0])
                e_ev.append(each_line[1])
                f.append(each_line[2])
                s.append(each_line[3])
        f_max = np.amax(f)
        df['E_eV'] = e_ev
        df['t_us'] = t_us
        df['f'] = f
        df['s'] = s
        df['f_norm'] = f/f_max
        df['s_norm'] = s/f_max
        # df[]
        # file_name = 'energy_' + str(index + 1) + '.csv'
        # df.to_csv(file_name, index=False)
        shape_dict[each_energy] = df
    return shape_dict
