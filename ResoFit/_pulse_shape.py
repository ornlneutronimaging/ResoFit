import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit import Model
from lmfit.models import LinearModel

import ResoFit._utilities as fit_util
from ResoFit.model import cole_windsor
from ResoFit.model import cole_windsor_jparc
from ResoFit.model import ikeda_carpenter
from ResoFit.model import ikeda_carpenter_jparc
from ResoFit.model import loglog_linear


class NeutronPulse(object):

    def __init__(self, path, model_index=1):
        """
        Load total neutron pulse shape from .dat file

        :param path: path to the '.dat' file
        :type path: str
        """
        self.shape_total_df = load_neutron_total_shape(path)
        self.shape_dict = None
        self.result_shape_fit = None
        # self.result_param_fit = None
        self.fitted_param_df_dir = None
        self.param_df = None
        self.model_param_names = None
        self.e_min = None
        self.e_max = None
        self.result_neutron_folder = None
        self.model_map = {1: 'ikeda_carpenter',
                          2: 'cole_windsor',
                          3: 'pseudo_voigt',
                          4: 'ikeda_carpenter_jparc',
                          5: 'cole_windsor_jparc'
                          }
        self.model_index = model_index
        self.model_used = self.model_map[model_index]
        if self.result_neutron_folder is None:
            self.result_neutron_folder = self._check_and_make_subdir('result', 'neutron_pulse', self.model_used)
        # self.model_params = None
        # self.params_to_fitshape = None
        # self.shape_result = None
        # self.model_index = None

    def plot_total(self, x_type='energy'):
        x_type_list = ['energy', 'lambda']
        if x_type not in x_type_list:
            raise ValueError("Please specify the x-axis type using one from '{}'.".format(x_type_list))
        if x_type == 'energy':
            plt.loglog(self.shape_total_df['E_eV'], self.shape_total_df['f(E)'], '.')
            plt.xlabel('Energy (eV)')
            plt.ylabel('Flux (n/sterad/pulse/eV)')
        elif x_type == 'lambda':
            plt.loglog(self.shape_total_df['I_angstrom'], self.shape_total_df['f(I)'], '.')
            plt.xlabel(u"Wavelength (\u212B)")
            plt.ylabel('Flux (n/sterad/pulse/Angstrom)')

    def load_shape_each(self, path, save_each=False):
        """
        Load each eV neutron pulse shape from .dat file

        :param save_each:
        :type save_each:
        :param path: path to the '.dat' file
        :type path: str
        """
        self.shape_dict = load_neutron_each_shape(path, export=save_each)

    def fit_shape(self, e_min, e_max,
                  drop=False, norm=True, show_init=True,
                  check_each=False, save_fig=False, overwrite_csv=False):
        # [1: 'ikeda_carpenter', 2: 'cole_windsor', 3: 'pseudo_voigt']
        self.e_min = e_min
        self.e_max = e_max

        # File name of fitted_param_df.csv
        assert self.model_used is not None
        _e_min = str(self.e_min) + 'eV_'
        _e_max = str(self.e_max) + 'eV_'
        _model_s = self.model_used + '.csv'
        _filename = 'Neutron_fitted_params_' + _e_min + _e_max + _model_s
        self.fitted_param_df_dir = os.path.join(self.result_neutron_folder, _filename)

        # File exists
        if os.path.isfile(self.fitted_param_df_dir):
            print("'{}' exists...".format(self.fitted_param_df_dir))
            if overwrite_csv:
                # Override==True, perform fitting and overwrite the .csv file
                print("File overwriting...")
                print("New fitting starts...")
                # Fitting starts
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
                        # f = self.shape_dict[each_e][_data_used]['f']
                        t = self.shape_dict[each_e][_data_used]['t_us']
                        param_dict_fitted[each_e]['fitted_params'] = self._fit_shape(f=f,
                                                                                     t=t,
                                                                                     e=each_e,
                                                                                     model_index=self.model_index,
                                                                                     check_each=check_each,
                                                                                     show_init=show_init,
                                                                                     save_fig=save_fig)
                print("File overwritten.")
                # Organize fitted parameters into pd.DataFrame
                self.param_df = self._form_params_df(param_dict=param_dict_fitted, save=True)
            else:
                # Override==False, read the .csv file
                self.param_df = pd.read_csv(self.fitted_param_df_dir)
                print("File loaded.")

        # File not exists, perform fitting
        else:
            print("No previous fitting file detected.\nNew fitting starts...")
            # Fitting starts
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
                    if norm:
                        _flux_used = 'f_norm'
                    else:
                        _flux_used = 'f'
                    f = self.shape_dict[each_e][_data_used][_flux_used]
                    t = self.shape_dict[each_e][_data_used]['t_us']
                    param_dict_fitted[each_e]['fitted_params'] = self._fit_shape(f=f,
                                                                                 t=t,
                                                                                 e=each_e,
                                                                                 model_index=self.model_index,
                                                                                 check_each=check_each,
                                                                                 show_init=show_init,
                                                                                 save_fig=save_fig)
            # Organize fitted parameters into pd.DataFrame
            self.param_df = self._form_params_df(param_dict=param_dict_fitted, save=True)

    def plot_params_vs_e(self, loglog=True):
        assert self.param_df is not None
        self.param_df.set_index('E_eV').plot(loglog=loglog, style='.')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Fitted parameter value')

    def fit_params(self, show_init=True, check_each=False, save_fig=False):
        if self.param_df is None:
            raise ValueError("'NeutronPulse.fit_shape' must be applied before 'NeutronPulse.fit_params'")
        _param_df = self.param_df
        # param_df_log10 = np.log10(_param_df)
        e_log = np.log10(_param_df['E_eV'])
        param_name_list = list(_param_df.columns.drop('E_eV'))
        param_dict_linear = {}
        for each_param in param_name_list:
            # param_dict_linear[each_param] = self._fit_params(y=param_df_log10[each_param],
            #                                                  x=e_log,
            #                                                  name=each_param,
            #                                                  check_each=check_each,
            #                                                  show_init=show_init,
            #                                                  save_fig=save_fig
            #                                                  )
            param_dict_linear[each_param] = self._fit_params(y=_param_df[each_param],
                                                             x=_param_df['E_eV'],
                                                             x_log=e_log,
                                                             name=each_param,
                                                             check_each=check_each,
                                                             show_init=show_init,
                                                             save_fig=save_fig
                                                             )

        return param_dict_linear

    def _fit_params(self, y, x, x_log, name, show_init, check_each, save_fig):

        # my_param_model = LinearModel()
        # assert my_param_model.independent_vars == ['x']
        #
        # # Guess params
        # params = my_param_model.guess(y, x)
        # # Make params
        #
        # # Fit the model
        # out = my_param_model.fit(y, params, x=x)

        # Guess  to make param
        y_log = np.log10(y)
        _temp_guess_model = LinearModel()
        params = _temp_guess_model.guess(y_log, x_log)

        # Load actual loglog linear model
        my_param_model = Model(loglog_linear)
        assert my_param_model.independent_vars == ['x']
        assert my_param_model.param_names == ['slope', 'intercept']

        # Fit the model
        out = my_param_model.fit(y, params, x=x)

        if check_each:
            # assert self.result_shape_fit is not None
            if name == 'f_max':
                _y_label = 'Flux (neutrons/s/cm2)'
            else:
                _y_label = 'Parameter value (arb. unit)'
            _x_label = 'E (eV)'
            if show_init is True:
                out.plot(xlabel=_x_label, ylabel=_y_label)
            else:
                out.plot(xlabel=_x_label, ylabel=_y_label, initfmt='None')

            plt.title(name, fontsize=12)
            plt.yscale('log')
            plt.xscale('log')
            plt.xlim(left=self.e_min * 0.8, right=self.e_max * 1.2)
            plt.ylim(bottom=np.amin(y) * 0.8, top=np.amax(y) * 1.2)

            if save_fig:
                # Check and make dir to save
                assert self.result_neutron_folder is not None
                assert self.model_used is not None
                _model_s = self.model_used + '_' + name + '.png'
                _filename = 'Neutron_pulse_fit_' + _model_s
                _dir_to_save = os.path.join(self.result_neutron_folder, _filename)
                plt.savefig(_dir_to_save, dpi=300, transparent=False)
                plt.close()
            else:
                plt.show()
        elif save_fig:
            raise ValueError("'check_each' has to be 'True' in order to save figure")

        return out.best_values

    def _fit_shape(self, f, t, e, model_index, show_init, check_each, save_fig):
        _model_map = self.model_map
        verbose = False

        # ikeda_carpenter
        if model_index == 1:
            self.model_used = _model_map[model_index]
            my_model = Model(ikeda_carpenter)

            assert my_model.param_names == ['alpha', 'beta', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']
            self.model_param_names = my_model.param_names

            # Set params hints
            if self.result_shape_fit is None:
                my_model.set_param_hint('alpha', value=0.699, min=0, max=20)
                my_model.set_param_hint('beta', value=0.0215, min=0, max=20)
                my_model.set_param_hint('fraction', value=0.383, min=0, max=1)
                my_model.set_param_hint('t0', value=0.0889, min=0, max=20)
                my_model.set_param_hint('norm_factor', value=1, min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)
            else:
                _best_pre_values = self.result_shape_fit.best_values
                my_model.set_param_hint('alpha', value=_best_pre_values['alpha'], min=0, max=20)
                my_model.set_param_hint('beta', value=_best_pre_values['beta'], min=0, max=20)
                my_model.set_param_hint('fraction', value=_best_pre_values['fraction'], min=0, max=1)
                my_model.set_param_hint('t0', value=_best_pre_values['t0'], min=0, max=20)
                my_model.set_param_hint('norm_factor', value=_best_pre_values['norm_factor'], min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)

            # Fit the model
            self.result_shape_fit = my_model.fit(f, params, t=t)

        # ikeda_carpenter_jparc
        elif model_index == 4:
            self.model_used = _model_map[model_index]
            my_model = Model(ikeda_carpenter_jparc)

            assert my_model.param_names == ['alpha', 'beta', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']
            self.model_param_names = my_model.param_names

            # Set params hints
            if self.result_shape_fit is None:
                my_model.set_param_hint('alpha', value=0.699, min=0, max=20)
                my_model.set_param_hint('beta', value=0.0215, min=0, max=20)
                my_model.set_param_hint('fraction', value=0.383, min=0, max=1)
                my_model.set_param_hint('t0', value=0.0889, min=0, max=20)
                my_model.set_param_hint('norm_factor', value=1, min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)
            else:
                _best_pre_values = self.result_shape_fit.best_values
                my_model.set_param_hint('alpha', value=_best_pre_values['alpha'], min=0, max=20)
                my_model.set_param_hint('beta', value=_best_pre_values['beta'], min=0, max=20)
                my_model.set_param_hint('fraction', value=_best_pre_values['fraction'], min=0, max=1)
                my_model.set_param_hint('t0', value=_best_pre_values['t0'], min=0, max=20)
                my_model.set_param_hint('norm_factor', value=_best_pre_values['norm_factor'], min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)

            # Fit the model
            self.result_shape_fit = my_model.fit(f, params, t=t)

        # cole_windsor
        elif model_index == 2:
            self.model_used = _model_map[model_index]
            my_model = Model(cole_windsor)

            assert my_model.param_names == ['sig1', 'sig2', 'gamma', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']
            self.model_param_names = my_model.param_names

            # Set params hints
            if self.result_shape_fit is None:
                my_model.set_param_hint('sig1', value=0.06917, min=0, max=20)
                my_model.set_param_hint('sig2', value=0.2041, min=0, max=20)
                my_model.set_param_hint('gamma', value=6.291, min=0, max=20)
                my_model.set_param_hint('fraction', value=0.1308, min=0, max=1)
                my_model.set_param_hint('t0', value=0.3176, min=0, max=20)
                my_model.set_param_hint('norm_factor', value=0.9951, min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)
            else:
                _best_pre_values = self.result_shape_fit.best_values
                my_model.set_param_hint('sig1', value=_best_pre_values['sig1'], min=0, max=20)
                my_model.set_param_hint('sig2', value=_best_pre_values['sig2'], min=0, max=20)
                my_model.set_param_hint('gamma', value=_best_pre_values['gamma'], min=0, max=20)
                my_model.set_param_hint('fraction', value=_best_pre_values['fraction'], min=0, max=1)
                my_model.set_param_hint('t0', value=_best_pre_values['t0'], min=0, max=20)
                my_model.set_param_hint('norm_factor', value=_best_pre_values['norm_factor'], min=0)

                # Make params
                params = my_model.make_params(verbose=verbose)

            # Fit the model
            self.result_shape_fit = my_model.fit(f, params, t=t)

        # cole_windsor_jparc
        elif model_index == 5:
            self.model_used = _model_map[model_index]
            my_model = Model(cole_windsor_jparc)

            assert my_model.param_names == ['sig1', 'sig2', 'gam1', 'gam2', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']
            self.model_param_names = my_model.param_names

            # Set params hints
            if self.result_shape_fit is None:
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
                _best_pre_values = self.result_shape_fit.best_values
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
            self.result_shape_fit = my_model.fit(f, params, t=t)

        else:
            raise ValueError("Model index not exists, please refer to: '{}' ".format(_model_map))

        # Check or visually inspect the fitting result of each energy
        if check_each:
            assert self.result_shape_fit is not None
            _y_label = 'Normalized neutron flux (arb. unit)'
            _x_label = 'Time (µs)'
            _e_text = str(e) + ' eV'
            if show_init is True:
                self.result_shape_fit.plot(xlabel=_x_label, ylabel=_y_label)
            else:
                self.result_shape_fit.plot(xlabel=_x_label, ylabel=_y_label, initfmt='None')
            plt.title(_e_text, fontsize=12)

            # set x_max for plot
            if e < 1:
                _plt_xmax = 40
            elif 1 <= e < 5:
                _plt_xmax = 20
            elif 5 <= e < 15:
                _plt_xmax = 10
            elif 15 <= e < 30:
                _plt_xmax = 6
            elif 30 <= e < 50:
                _plt_xmax = 4
            elif 50 <= e < 500:
                _plt_xmax = 2.5
            else:
                _plt_xmax = 1.5

            # Show fitted params in figure
            _param_text = ''
            for _each in self.model_param_names:
                _param_text = _param_text + _each + ': ' + str(self.result_shape_fit.best_values[_each]) + '\n'
            param_text = _param_text[:-2]  # get rid of '\n' at the end
            plt.text(x=0.5 * _plt_xmax, y=0.4, s=param_text, fontsize=10,
                     bbox={'facecolor': 'None', 'alpha': 0.5, 'pad': 2}
                     )

            # Apply plot limit
            plt.xlim(xmin=0, xmax=_plt_xmax)

            if save_fig:
                # Check and make dir to save
                assert self.result_neutron_folder is not None
                _model_s = _model_map[model_index] + '.png'
                _filename = 'Neutron_pulse_fit_' + str(e) + '_eV_' + _model_s
                _dir_to_save = os.path.join(self.result_neutron_folder, _filename)
                plt.savefig(_dir_to_save, dpi=300, transparent=False)
                plt.close()
            else:
                plt.show()
        elif save_fig:
            raise ValueError("'check_each' has to be 'True' in order to save figure")

        return self.result_shape_fit.best_values

    def _form_params_df(self, param_dict, save=False):
        e_list_fitted = list(param_dict.keys())
        if self.model_param_names is None:
            param_list_fitted = list(param_dict[e_list_fitted[0]]['fitted_params'].keys())
        else:
            param_list_fitted = self.model_param_names

        _dict = {}
        for each_param in param_list_fitted:
            _dict[each_param] = []
        f_max_list = []
        for each_e in e_list_fitted:
            f_max_list.append(self.shape_dict[each_e]['f_max'])
            for each_param in param_list_fitted:
                _dict[each_param].append(param_dict[each_e]['fitted_params'][each_param])
        _df = pd.DataFrame()
        _df['E_eV'] = e_list_fitted
        for each_param in param_list_fitted:
            _df[each_param] = _dict[each_param]
        _df['f_max'] = f_max_list
        print(_df)
        print("NOTE: 'f_max' in the params_df is NOT a fitted parameter, it is the maximum in each shape.")

        if save is True:
            assert self.result_neutron_folder is not None
            assert self.fitted_param_df_dir is not None
            _df.to_csv(path_or_buf=self.fitted_param_df_dir, index=False)

        return _df

    def _export_total(self, filename=None):
        assert self.shape_total_df is not None

        if filename is None:
            self.shape_total_df.to_clipboard(excel=True)
        else:
            self.shape_total_df.to_csv(filename)

    def _export_each_to_csv(self):
        for index, each_energy in enumerate(self.shape_dict.keys()):
            df = self.shape_dict[each_energy]
            file_name = 'energy_' + str(index + 1) + '.csv'
            df.to_csv(file_name, index=False)
            print("Neutron pulse shape of 'E = {} eV' has exported to './{}'".format(each_energy, file_name))

    @staticmethod
    def _check_and_make_subdir(*args):
        # Check and make dir to save
        _file_path = os.path.abspath(os.path.dirname(__file__))
        for arg in args:
            _created_path = fit_util.check_and_make_dir(_file_path, arg)
            _file_path = _created_path
        return _file_path


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


def load_neutron_each_shape(path, export=False):
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
        df['E_eV'] = e_ev
        data_dict['data_raw'] = df
        df_dropped = df.drop(df[df.f_norm <= 0].index)
        df_dropped.reset_index(drop=True, inplace=True)
        data_dict['data'] = df_dropped
        data_dict['f_max'] = f_max
        shape_dict[each_energy] = data_dict

        # Export to .csv
        if export is True:
            file_name = 'energy_' + str(index + 1) + '.csv'
            df.to_csv(file_name, index=False)

    return shape_dict


def load_proton_pulse(path):
    df = pd.read_csv(path, sep=' ', skiprows=1, header=None)
    df.columns = ['t_ns', 'intensity']
    return df
