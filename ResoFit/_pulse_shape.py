import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import Parameters
from lmfit import minimize
from lmfit import Model
from ResoFit.model import cole_windsor
from ResoFit.model import cole_windsor_jparc
from ResoFit.model import ikeda_carpenter
import ImagingReso._utilities as reso_util
import ResoFit._utilities as fit_util


class NeutronPulse(object):

    def __init__(self, path):
        """"""
        self.shape_total_df = load_neutron_total_shape(path)
        # self.params_to_fitshape = None
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

    def export_each_to_csv(self):
        for index, each_energy in enumerate(self.shape_dict.keys()):
            df = self.shape_dict[each_energy]
            file_name = 'energy_' + str(index + 1) + '.csv'
            df.to_csv(file_name, index=False)
            print("Neutron pulse shape of 'E = {} eV' has exported to './{}'".format(each_energy, file_name))

    def fit_shape(self, model_index):
        # [1: 'ikeda_carpenter', 2: 'cole_windsor', 3: 'pseudo_voigt']
        for each in self.shape_dict.keys():
            f = self.shape_dict[each]['data']['f_norm']
            t = self.shape_dict[each]['data']['t_us']
            self.shape_dict[each]['fitted_params'] = self._fit_shape(f=f, t=t, model_index=model_index)

    def fit_params(self):
        pass

    def _form_fitted_df(self):
        e_list = list(self.shape_dict.keys())
        assert 'fitted_params' in self.shape_dict[e_list[0]]
        params_list = list(self.shape_dict[e_list[0]]['fitted_params'].keys())
        _dict = {}
        for each_param in params_list:
            _dict[each_param] = []
        for each_e in e_list:
            for each_param in params_list:
                _dict[each_param].append(self.shape_dict[each_e]['fitted_params'][each_param])
        _df = pd.DataFrame()
        for each_param in params_list:
            _df[each_param] = _dict[each_param]

        return _df

    def _fit_shape(self, f, t, model_index):
        # [1: 'ikeda_carpenter', 2: 'cole_windsor', 3: 'pseudo_voigt']
        self.model_index = model_index
        # self.params_to_fitshape = Parameters()

        # ikeda_carpenter
        if self.model_index == 1:
            my_model = Model(ikeda_carpenter)
            model_params = ['alpha', 'beta', 'fraction', 't0', 'norm_factor']
            assert my_model.param_names == ['alpha', 'beta', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']

            # Set params hints
            my_model.set_param_hint('alpha', value=0.699, min=0, max=20)
            my_model.set_param_hint('beta', value=0.0215, min=0, max=20)
            my_model.set_param_hint('fraction', value=0.383, min=0, max=1)
            my_model.set_param_hint('t0', value=0.0889, min=0, max=20)
            my_model.set_param_hint('norm_factor', value=1, min=0)

            # Make params
            params = my_model.make_params(verbose=True)

            # Fit the model
            result = my_model.fit(f, params, t=t)

            # result.best_values

        # cole_windsor
        elif self.model_index == 2:
            my_model = Model(cole_windsor)

            assert my_model.param_names == ['sig1', 'sig2', 'gamma', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']

            # Set params hints
            my_model.set_param_hint('sig1', value=0.06917, min=0, max=20)
            my_model.set_param_hint('sig2', value=0.2041, min=0, max=20)
            my_model.set_param_hint('gamma', value=6.291, min=0, max=20)
            my_model.set_param_hint('fraction', value=0.1308, min=0, max=1)
            my_model.set_param_hint('t0', value=0.3176, min=0, max=20)
            my_model.set_param_hint('norm_factor', value=0.9951, min=0)

            # Make params
            params = my_model.make_params(verbose=True)

            # Fit the model
            result = my_model.fit(f, params, t=t)

        # cole_windsor_jparc
        elif self.model_index == 3:
            my_model = Model(cole_windsor_jparc)

            assert my_model.param_names == ['sig1', 'sig2', 'gam1', 'gam2', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']

            # Set params hints
            my_model.set_param_hint('sig1', value=0.06917, min=0, max=20)
            my_model.set_param_hint('sig2', value=0.2041, min=0, max=20)
            my_model.set_param_hint('gam1', value=6.291, min=0, max=20)
            my_model.set_param_hint('gam2', value=1.285, min=0, max=20)
            my_model.set_param_hint('fraction', value=0.1308, min=0, max=1)
            my_model.set_param_hint('t0', value=0.3176, min=0, max=20)
            my_model.set_param_hint('norm_factor', value=0.9951, min=0)

            # Make params
            params = my_model.make_params(verbose=True)

            # Fit the model
            result = my_model.fit(f, params, t=t)

        return result.best_values

        # # Print before
        # print("+----------------- Fit neutron pulse shape -----------------+\nParams before:")
        # self.params_to_fitshape.pretty_print()
        # # Use lmfit to obtain params by minimizing gap_function
        #
        # # Print after
        # print("\nParams after:")
        # self.shape_result.__dict__['params'].pretty_print()
        # # Print chi^2
        # print("Calibration chi^2 : {}\n".format(self.shape_result.__dict__['chisqr']))


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
        data_dict = {}
        t_us = []
        # e_ev = []
        f = []
        s = []
        df = pd.DataFrame()
        for each_line in q:
            if each_energy == each_line[1]:
                t_us.append(each_line[0])
                # e_ev.append(each_line[1])
                f.append(each_line[2])
                s.append(each_line[3])
        f_max = np.amax(f)
        # df['E_eV'] = e_ev
        df['t_us'] = t_us
        df['f'] = f
        df['s'] = s
        df['f_norm'] = f / f_max
        df['s_norm'] = s / f_max
        # df[]
        # file_name = 'energy_' + str(index + 1) + '.csv'
        # df.to_csv(file_name, index=False)
        data_dict['data'] = df
        shape_dict[each_energy] = data_dict
    return shape_dict
