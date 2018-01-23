import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import Parameters
from lmfit import minimize
from lmfit import Model
from ResoFit.model import cole_windsor
from ResoFit.model import cole_windsor_jparc
from ResoFit.model import ikeda_carpenter
from ResoFit.model import ikeda_carpenter_jparc
import ImagingReso._utilities as reso_util
import ResoFit._utilities as fit_util
import os


class NeutronPulse(object):

    def __init__(self, path):
        """"""
        self.shape_total_df = load_neutron_total_shape(path)
        self.shape_dict = None
        self.result = None
        self.param_df_fitted = None
        self.model_used = None
        self.model_param_names = None
        self.e_min = None
        self.e_max = None
        self.result_neutron_folder = None
        # self.model_params = None
        # self.params_to_fitshape = None
        # self.shape_result = None
        # self.model_index = None

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

    def fit_shape(self, e_min, e_max, model_index=1, drop=False, show_init=True, check_each=False, save_fig=False, save_df=False):
        # [1: 'ikeda_carpenter', 2: 'cole_windsor', 3: 'pseudo_voigt']
        self.e_min = e_min
        self.e_max = e_max
        param_dict_fitted = {}
        for each_e in self.shape_dict.keys():
            if e_min <= each_e <= e_max:
                print("Fitting [{} eV] ...".format(each_e))
                param_dict_fitted[each_e] = {}
                # If drop is True, only flux value above 0 will be used
                if drop:
                    _data_used = 'data'
                else:
                    _data_used = 'data_raw'
                f = self.shape_dict[each_e][_data_used]['f_norm']
                t = self.shape_dict[each_e][_data_used]['t_us']
                param_dict_fitted[each_e]['fitted_params'] = self._fit_shape(f=f,
                                                                             t=t,
                                                                             e=each_e,
                                                                             model_index=model_index,
                                                                             check_each=check_each,
                                                                             show_init=show_init,
                                                                             save_fig=save_fig
                                                                             )
        # Organize fitted parameters into pd.DataFrame
        self.param_df_fitted = self._form_fitted_df(param_dict=param_dict_fitted, save=save_df)
        print(self.param_df_fitted)

    def fit_params(self):
        if self.param_df_fitted is None:
            raise ValueError("'.fit_shape' must be applied before '.fit_params'")

        pass

    def _form_fitted_df(self, param_dict, save=False):
        e_list_fitted = list(param_dict.keys())
        if self.model_param_names is None:
            param_list_fitted = list(param_dict[e_list_fitted[0]]['fitted_params'].keys())
        else:
            param_list_fitted = self.model_param_names

        _dict = {}
        for each_param in param_list_fitted:
            _dict[each_param] = []
        for each_e in e_list_fitted:
            for each_param in param_list_fitted:
                _dict[each_param].append(param_dict[each_e]['fitted_params'][each_param])
        _df = pd.DataFrame()
        _df['E_eV'] = e_list_fitted
        for each_param in param_list_fitted:
            _df[each_param] = _dict[each_param]

        if save is True:
            # Check and make dir to save
            if self.result_neutron_folder is None:
                self.result_neutron_folder = self._check_and_make_subdir('result', 'neutron_pulse', self.model_used)

            # Form file name
            _e_min = str(self.e_min) + 'eV_'
            _e_max = str(self.e_max) + 'eV_'
            _model_s = self.model_used + '.csv'
            _file_name = 'Neutron_fitted_params_' + _e_min + _e_max + _model_s

            # Save
            _dir_to_save = os.path.join(self.result_neutron_folder, _file_name)
            _df.to_csv(path_or_buf=_dir_to_save, index=False)

        return _df

    @staticmethod
    def _check_and_make_subdir(*args):
        # Check and make dir to save
        _file_path = os.path.abspath(os.path.dirname(__file__))
        for arg in args:
            _created_path = fit_util.check_and_make_dir(_file_path, arg)
            _file_path = _created_path
        return _file_path

    def _fit_shape(self, f, t, e, model_index, show_init, check_each, save_fig):
        _model_map = {1: 'ikeda_carpenter',
                      2: 'cole_windsor',
                      3: 'pseudo_voigt',
                      4: 'ikeda_carpenter_jparc',
                      5: 'cole_windsor_jparc',
                      }
        # [1: 'ikeda_carpenter', 2: 'cole_windsor', 3: 'pseudo_voigt']
        verbose = False

        # ikeda_carpenter
        if model_index == 1:
            self.model_used = _model_map[model_index]
            my_model = Model(ikeda_carpenter)

            assert my_model.param_names == ['alpha', 'beta', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']
            self.model_param_names = my_model.param_names

            # Set params hints
            if self.result is None:
                my_model.set_param_hint('alpha', value=0.699, min=0, max=20)
                my_model.set_param_hint('beta', value=0.0215, min=0, max=20)
                my_model.set_param_hint('fraction', value=0.383, min=0, max=1)
                my_model.set_param_hint('t0', value=0.0889, min=0, max=20)
                my_model.set_param_hint('norm_factor', value=1, min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)
            else:
                _best_pre_values = self.result.best_values
                my_model.set_param_hint('alpha', value=_best_pre_values['alpha'], min=0, max=20)
                my_model.set_param_hint('beta', value=_best_pre_values['beta'], min=0, max=20)
                my_model.set_param_hint('fraction', value=_best_pre_values['fraction'], min=0, max=1)
                my_model.set_param_hint('t0', value=_best_pre_values['t0'], min=0, max=20)
                my_model.set_param_hint('norm_factor', value=_best_pre_values['norm_factor'], min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)

            # Fit the model
            self.result = my_model.fit(f, params, t=t)

        # ikeda_carpenter_jparc
        elif model_index == 4:
            self.model_used = _model_map[model_index]
            my_model = Model(ikeda_carpenter_jparc)

            assert my_model.param_names == ['alpha', 'beta', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']
            self.model_param_names = my_model.param_names

            # Set params hints
            if self.result is None:
                my_model.set_param_hint('alpha', value=0.699, min=0, max=20)
                my_model.set_param_hint('beta', value=0.0215, min=0, max=20)
                my_model.set_param_hint('fraction', value=0.383, min=0, max=1)
                my_model.set_param_hint('t0', value=0.0889, min=0, max=20)
                my_model.set_param_hint('norm_factor', value=1, min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)
            else:
                _best_pre_values = self.result.best_values
                my_model.set_param_hint('alpha', value=_best_pre_values['alpha'], min=0, max=20)
                my_model.set_param_hint('beta', value=_best_pre_values['beta'], min=0, max=20)
                my_model.set_param_hint('fraction', value=_best_pre_values['fraction'], min=0, max=1)
                my_model.set_param_hint('t0', value=_best_pre_values['t0'], min=0, max=20)
                my_model.set_param_hint('norm_factor', value=_best_pre_values['norm_factor'], min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)

            # Fit the model
            self.result = my_model.fit(f, params, t=t)

        # cole_windsor
        elif model_index == 2:
            self.model_used = _model_map[model_index]
            my_model = Model(cole_windsor)

            assert my_model.param_names == ['sig1', 'sig2', 'gamma', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']
            self.model_param_names = my_model.param_names

            # Set params hints
            if self.result is None:
                my_model.set_param_hint('sig1', value=0.06917, min=0, max=20)
                my_model.set_param_hint('sig2', value=0.2041, min=0, max=20)
                my_model.set_param_hint('gamma', value=6.291, min=0, max=20)
                my_model.set_param_hint('fraction', value=0.1308, min=0, max=1)
                my_model.set_param_hint('t0', value=0.3176, min=0, max=20)
                my_model.set_param_hint('norm_factor', value=0.9951, min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)
            else:
                _best_pre_values = self.result.best_values
                my_model.set_param_hint('sig1', value=_best_pre_values['sig1'], min=0, max=20)
                my_model.set_param_hint('sig2', value=_best_pre_values['sig2'], min=0, max=20)
                my_model.set_param_hint('gamma', value=_best_pre_values['gamma'], min=0, max=20)
                my_model.set_param_hint('fraction', value=_best_pre_values['fraction'], min=0, max=1)
                my_model.set_param_hint('t0', value=_best_pre_values['t0'], min=0, max=20)
                my_model.set_param_hint('norm_factor', value=_best_pre_values['norm_factor'], min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)

            # Fit the model
            self.result = my_model.fit(f, params, t=t)

        # cole_windsor_jparc
        elif model_index == 5:
            self.model_used = _model_map[model_index]
            my_model = Model(cole_windsor_jparc)

            assert my_model.param_names == ['sig1', 'sig2', 'gam1', 'gam2', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']
            self.model_param_names = my_model.param_names

            # Set params hints
            if self.result is None:
                my_model.set_param_hint('sig1', value=0.06917, min=0, max=20)
                my_model.set_param_hint('sig2', value=0.2041, min=0, max=20)
                my_model.set_param_hint('gam1', value=6.291, min=0, max=20)
                my_model.set_param_hint('gam2', value=1.285, min=0, max=20)
                my_model.set_param_hint('fraction', value=0.1308, min=0, max=1)
                my_model.set_param_hint('t0', value=0.3176, min=0, max=20)
                my_model.set_param_hint('norm_factor', value=0.9951, min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)
            else:
                _best_pre_values = self.result.best_values
                my_model.set_param_hint('sig1', value=_best_pre_values['sig1'], min=0, max=20)
                my_model.set_param_hint('sig2', value=_best_pre_values['sig2'], min=0, max=20)
                my_model.set_param_hint('gam1', value=_best_pre_values['gam1'], min=0, max=20)
                my_model.set_param_hint('gam2', value=_best_pre_values['gam2'], min=0, max=20)
                my_model.set_param_hint('fraction', value=_best_pre_values['fraction'], min=0, max=1)
                my_model.set_param_hint('t0', value=_best_pre_values['t0'], min=0, max=20)
                my_model.set_param_hint('norm_factor', value=_best_pre_values['norm_factor'], min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)

            # Fit the model
            self.result = my_model.fit(f, params, t=t)

        else:
            raise ValueError("Model index not exists, please refer to: '{}' ".format(_model_map))

        if check_each:
            assert self.result is not None
            _y_label = 'Normalized neutron flux (arb. unit)'
            _x_label = 'Time (µs)'
            _e_text = str(e) + ' eV'
            if show_init is True:
                self.result.plot(xlabel=_x_label, ylabel=_y_label)
            else:
                self.result.plot(xlabel=_x_label, ylabel=_y_label, initfmt='None')
            plt.text(x=0, y=1.09, s=_e_text, fontsize=12)

            if e < 1:
                plt.xlim(xmin=0, xmax=40)
            elif 1 <= e < 5:
                plt.xlim(xmin=0, xmax=20)
            elif 5 <= e < 15:
                plt.xlim(xmin=0, xmax=10)
            elif 15 <= e < 30:
                plt.xlim(xmin=0, xmax=6)
            elif 30 <= e < 50:
                plt.xlim(xmin=0, xmax=4)
            elif 50 <= e < 500:
                plt.xlim(xmin=0, xmax=2.5)
            else:
                plt.xlim(xmin=0, xmax=1.5)

            if save_fig:
                # Check and make dir to save
                if self.result_neutron_folder is None:
                    self.result_neutron_folder = self._check_and_make_subdir('result', 'neutron_pulse', self.model_used)

                _model_s = _model_map[model_index] + '.png'
                _filename = 'Neutron_pulse_fit_' + str(e) + '_eV_' + _model_s
                _dir_to_save = os.path.join(self.result_neutron_folder, _filename)
                plt.savefig(_dir_to_save, dpi=300, transparent=False)
                plt.close()
            else:
                plt.show()
        elif save_fig:
            raise ValueError("'check_each' has to be 'True' in order to save figure")

        return self.result.best_values

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

    def __init__(self, path):
        """"""
        self.shape_df = load_proton_pulse(path)

    def fit_shape(self):
        t_ns = self.shape_df['t_ns']
        intensity = self.shape_df['intensity']
        # my_model = Model(guass)
        # result = my_model.fit((intensity, params, t=t_ns))
        pass


# Functions to load files #
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
        df['t_us'] = t_us
        df['f'] = f
        df['s'] = s
        df['f_norm'] = f / f_max
        df['s_norm'] = s / f_max
        # df['f_max'] = [f_max] * len(f)
        df['E_eV'] = e_ev
        data_dict['data_raw'] = df
        df_dropped = df.drop(df[df.f_norm <= 0].index)
        df_dropped.reset_index(drop=True, inplace=True)
        data_dict['data'] = df_dropped
        shape_dict[each_energy] = data_dict

        # Export to .csv
        # df[]
        # file_name = 'energy_' + str(index + 1) + '.csv'
        # df.to_csv(file_name, index=False)
    return shape_dict


def load_proton_pulse(path):
    df = pd.read_csv(path, sep=' ', skiprows=1, header=None)
    df.columns = ['t_ns', 'intensity']
    return df
