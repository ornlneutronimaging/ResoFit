import os

import lmfit
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ImagingReso._utilities import ev_to_s
from cycler import cycler
from lmfit import Model
from lmfit.models import LinearModel
from scipy.interpolate import interp1d

import ResoFit._utilities as fit_util
from ResoFit.model import cole_windsor
from ResoFit.model import cole_windsor_jparc
from ResoFit.model import ikeda_carpenter
from ResoFit.model import ikeda_carpenter_jparc
from ResoFit.model import loglog_linear
from ResoFit.model import pseudo_voigt

# import ImagingReso._utilities as reso_util

__file_path = os.path.abspath(os.path.dirname(__file__))
_rel_path_to_proton = 'data/_data_for_tutorial/proton_pulse/waveform_20170901.txt'
proton_path = os.path.join(__file_path, _rel_path_to_proton)

t_min_us = 5e-2
t_max_us = 183.7
t_step_us = 0.01
t_nbr = round((t_max_us - t_min_us) / t_step_us) + 1

mcnp_plot_lim_dict = {'x_min': 0,
                      'x_max': 5,
                      'y_min': -0.05,
                      'y_max': 1.05}
plot_x_range = mcnp_plot_lim_dict['x_max'] - mcnp_plot_lim_dict['x_min']
default_cycler = cycler('color',
                        ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                         '#bcbd22', '#17becf'])
# fmt_list = ['o--', 'v--', 's--', 'p--', 'D--', 'h--', '^--', '<--', '>--', 'o--']
marker_list = ['o', 'v', 's', '^', 'D', 'h', '<', '>', 'x', 'p']


class NeutronPulse(object):

    def __init__(self, path, model_index=1):
        """
        Load total neutron pulse shape from .dat file

        :param path: path to the '.dat' file
        :type path: str
        """
        self.shape_total_df = _load_neutron_total_shape(path)
        self.shape_dict_mcnp = None
        self.shape_df_mcnp = None
        self.shape_df_mcnp_norm = None
        self.shape_df_interp = None
        self.shape_tof_df_interp = None
        self.shape_tof_df_dir = None

        self.result_shape_fit = None
        self.param_df_dir = None
        self.param_df = None
        self.linear_df = None
        self.linear_df_dir = None
        self.model_param_names = None
        self.e_min = None
        self.e_max = None
        self.t_us_mcnp = None
        self.t_us_conv_proton = np.linspace(t_min_us, t_max_us, t_nbr).round(3)
        self.result_neutron_folder = None
        self._energy_list = None
        self._energy_list_dropped = None

        self.model_map = {1: 'ikeda_carpenter',
                          2: 'cole_windsor',
                          3: 'pseudo_voigt',
                          4: 'ikeda_carpenter_jparc',
                          5: 'cole_windsor_jparc',
                          }
        self.model_index = None
        self.model_used = None
        self.model = None
        self.__set_model(model_index)

        if self.result_neutron_folder is None:
            self.result_neutron_folder = self._check_and_make_subdir('result', 'neutron_pulse', self.model_used)

        self.proton_pulse = ProtonPulse(path=proton_path)

    def load_shape_each(self, path, save_each=False):
        """
        Load each eV neutron pulse shape from .dat file

        :param save_each:
        :type save_each:
        :param path: path to the '.dat' file
        :type path: str
        """
        self.shape_dict_mcnp = _load_neutron_each_shape(path, export=save_each)
        self.shape_df_mcnp, self.shape_df_mcnp_norm = _shape_dict_to_dfs(self.shape_dict_mcnp,
                                                                         t_max=t_max_us)
        self._energy_list = list(self.shape_df_mcnp.set_index('t_us').columns)
        self.t_us_mcnp = np.array(self.shape_df_mcnp['t_us']).round(5)

    def _form_energy_list_dropped(self, e_min, e_max):
        assert self._energy_list is not None
        _energy_list_dropped = []
        for _each_e in self._energy_list:
            if e_min <= _each_e <= e_max:
                _energy_list_dropped.append(_each_e)
        _energy_list_dropped.sort()
        self._energy_list_dropped = _energy_list_dropped
        return _energy_list_dropped

    def plot_shape_total(self, x1_type='energy', x2_type='lambda', source_to_detector_m=None):
        """
        Plot the total beam shape obtained from MCNPX simulation

        :param x1_type:
        :type x1_type:
        :param x2_type:
        :type x2_type:
        :param source_to_detector_m:
        :type source_to_detector_m:
        :return: plot
        :rtype: matplotlib
        """
        x_type_list = ['energy', 'lambda', 'time', 'none']
        if x1_type not in x_type_list:
            raise ValueError("Please specify the x-axis type using one from '{}'.".format(x_type_list))
        if x2_type not in x_type_list:
            raise ValueError("Please specify the x-axis type using one from '{}'.".format(x_type_list))
        if x1_type == 'time' or x2_type == 'time':
            if source_to_detector_m is None:
                raise ValueError("Please specify the source-to-detector distance in m.")
        if x1_type == x2_type:
            x2_type = 'none'

        fig, ax1 = plt.subplots()
        if x1_type == 'energy':
            ax1.loglog(self.shape_total_df['E_eV'], self.shape_total_df['f(E)'], 'b.')
            ax1.set_xlabel('Energy (eV)', color='b')
        elif x1_type == 'lambda':
            ax1.loglog(self.shape_total_df['l_angstrom'], self.shape_total_df['f(l)'], 'b.')
            ax1.set_xlabel(u"Wavelength (\u212B)", color='b')
            ax1.invert_xaxis()
        elif x1_type == 'time':
            ax1.loglog(ev_to_s(array=self.shape_total_df['E_eV'],
                               offset_us=0,
                               source_to_detector_m=source_to_detector_m) * 1e6,
                       self.shape_total_df['f(l)'], 'b.')
            ax1.set_xlabel(u"Time-of-flight (\u03BCs)", color='b')
            ax1.invert_xaxis()
        ax1.set_ylabel('Flux (n/sterad/pulse)')
        ax1.tick_params('x', colors='b', which='both')
        ax1.grid(axis='x', which='both', color='b', alpha=0.3)
        ax1.grid(axis='y', which='major', alpha=0.3)
        if x2_type != 'none':
            ax2 = ax1.twiny()
            if x2_type == 'energy':
                ax2.loglog(self.shape_total_df['E_eV'], self.shape_total_df['f(E)'], 'rx')
                ax2.set_xlabel('Energy (eV)', color='r')
            elif x2_type == 'lambda':
                ax2.loglog(self.shape_total_df['l_angstrom'], self.shape_total_df['f(l)'], 'rx')
                ax2.set_xlabel(u"Wavelength (\u212B)", color='r')
                ax2.invert_xaxis()
            elif x2_type == 'time':
                ax2.loglog(ev_to_s(array=self.shape_total_df['E_eV'],
                                   offset_us=0,
                                   source_to_detector_m=source_to_detector_m) * 1e6,
                           self.shape_total_df['f(l)'], 'rx')
                ax2.set_xlabel(u"Time-of-flight (\u03BCs)", color='r')
                ax2.invert_xaxis()
            ax2.grid(axis='x', which='both', color='r', alpha=0.3)
            ax2.tick_params('x', colors='r', which='both')
        # ax1.set_title('Neutron total flux', y=1.08, loc='left')
        return fig

    def plot_shape_mcnp(self, e_min, e_max, logy=False, norm=False, marker='o', ax_mpl=None, plt_arrow=True):
        """
        Plot each eV beam shape obtained from MCNPX simulation

        :param e_min:
        :type e_min:
        :param e_max:
        :type e_max:
        :param logy:
        :type logy:
        :param norm:
        :type norm:
        :param marker:
        :type marker:
        :param ax_mpl:
        :type ax_mpl:
        :param plt_arrow:
        :type plt_arrow:
        :return:
        :rtype:
        """
        assert self.shape_dict_mcnp is not None
        if norm:
            _shape_df = self.shape_df_mcnp_norm
            _y_label = 'Ratio out of max flux of each energy'
        else:
            _shape_df = self.shape_df_mcnp
            _y_label = 'Flux (n/sterad/pulse)'

        # Only the energy provided by MCNPX simulation will be filtered
        _energy_list_dropped = self._form_energy_list_dropped(e_min=e_min, e_max=e_max)
        if ax_mpl is None:
            fig, ax_mpl = plt.subplots()
        ax_mpl.set_prop_cycle(default_cycler)
        for each in _energy_list_dropped:
            if logy:
                ax_mpl.semilogy(_shape_df['t_us'],
                                _shape_df[each],
                                linestyle='-',
                                marker=marker,
                                # fmt=marker,
                                label=str(each) + ' eV (MCNPX)')
            else:
                ax_mpl.plot(_shape_df['t_us'],
                            _shape_df[each],
                            linestyle='-',
                            marker=marker,
                            # fmt=marker,
                            label=str(each) + ' eV (MCNPX)')
        ax_mpl.set_ylabel(_y_label)
        ax_mpl.set_xlabel(u'Time (\u03BCs)')
        ax_mpl.grid()
        if len(_energy_list_dropped) <= 7:
            ax_mpl.legend()
        else:
            if plt_arrow:
                _plot_ev_arrow_as_legend(ax=ax_mpl, ev_list=_energy_list_dropped)
        ax_mpl.set_title('Energy dependent neutron pulse shape (MCNPX)')
        ax_mpl.set_xlim(left=mcnp_plot_lim_dict['x_min'], right=mcnp_plot_lim_dict['x_max'])
        return ax_mpl

    def plot_shape_interp(self, e_ev, source_to_detector_m, conv_proton, proton_params={},
                          t_interp=None, logy=False, norm=False, for_sum=False,
                          marker='o', ax_mpl=None, plt_arrow=True):
        """
        Plot each eV beam shape obtained from the fitting approach

        :param e_ev:
        :type e_ev:
        :param source_to_detector_m:
        :type source_to_detector_m:
        :param conv_proton:
        :type conv_proton:
        :param proton_params:
        :type proton_params:
        :param t_interp:
        :type t_interp:
        :param logy:
        :type logy:
        :param norm:
        :type norm:
        :param for_sum:
        :type for_sum:
        :param marker:
        :type marker:
        :param ax_mpl:
        :type ax_mpl:
        :param plt_arrow:
        :type plt_arrow:
        :return:
        :rtype:
        """
        if t_interp is None:
            t_interp = self.t_us_mcnp
        self._make_shape(e_ev=e_ev, t_interp=t_interp, norm=norm, for_sum=for_sum,
                         source_to_detector_m=source_to_detector_m, print_tof=False,
                         conv_proton=conv_proton, proton_params=proton_params)
        _shape_df_interp = self.shape_df_interp

        _y_label = 'Flux (n/sterad/pulse)'
        if norm:
            _y_label = 'Ratio out of max flux of each energy'
        _for_sum_s = ''
        if for_sum:
            _for_sum_s = ' for_sum'
        _norm_s = ''
        if norm:
            _norm_s = ' norm'
        _conv_proton_s = ''
        if conv_proton:
            _conv_proton_s = ' proton'
        _proton_param_s = ''
        for _param in proton_params.keys():
            _proton_param_s = _proton_param_s + '_' + _param + '_' + str(proton_params[_param])
        _details_s = _norm_s + _conv_proton_s + _proton_param_s + _for_sum_s
        _title_s = 'Energy dependent neutron pulse shape (interp.{})'.format(_details_s)

        _energy_interp_list = list(_shape_df_interp.set_index('t_us').columns)
        _energy_interp_list.sort()
        if ax_mpl is None:
            fig, ax_mpl = plt.subplots()
        ax_mpl.set_prop_cycle(default_cycler)
        for each in _energy_interp_list:
            if logy:
                ax_mpl.semilogy(_shape_df_interp['t_us'],
                                _shape_df_interp[each],
                                linestyle='--',
                                marker=marker,
                                # fmt=fmt,
                                mfc='none',
                                label=str(each) + ' eV (interp.)')
            else:
                ax_mpl.plot(_shape_df_interp['t_us'],
                            _shape_df_interp[each],
                            linestyle='--',
                            marker=marker,
                            # fmt=fmt,
                            mfc='none',
                            label=str(each) + ' eV (interp.)')
        ax_mpl.set_ylabel(_y_label)
        ax_mpl.set_xlabel(u'Time (\u03BCs)')
        ax_mpl.grid()
        if len(_energy_interp_list) <= 7:
            ax_mpl.legend()
        else:
            if plt_arrow:
                _plot_ev_arrow_as_legend(ax=ax_mpl, ev_list=_energy_interp_list)

        ax_mpl.set_title(_title_s)
        return ax_mpl

    def plot_shape_each_compare(self, e_min, e_max, source_to_detector_m, conv_proton, proton_params={},
                                t_interp=None, logy=False, norm=False, for_sum=False):
        """
        Plot each eV beam shape obtained from MCNPX simulation and current fitting approach to compare

        :param e_min:
        :type e_min:
        :param e_max:
        :type e_max:
        :param source_to_detector_m:
        :type source_to_detector_m:
        :param conv_proton:
        :type conv_proton:
        :param proton_params:
        :type proton_params:
        :param t_interp:
        :type t_interp:
        :param logy:
        :type logy:
        :param norm:
        :type norm:
        :param for_sum:
        :type for_sum:
        :return:
        :rtype:
        """
        if t_interp is None:
            t_interp = self.t_us_mcnp
        _energy_list_dropped = self._form_energy_list_dropped(e_min=e_min, e_max=e_max)
        ax_mcnp = self.plot_shape_mcnp(e_min=e_min, e_max=e_max, norm=norm, logy=logy, plt_arrow=False)
        ax = self.plot_shape_interp(e_ev=_energy_list_dropped,
                                    source_to_detector_m=source_to_detector_m,
                                    conv_proton=conv_proton, proton_params=proton_params,
                                    t_interp=t_interp, logy=logy, norm=norm, for_sum=for_sum, ax_mpl=ax_mcnp)
        ax.set_title('Energy dependent neutron pulse shape MCNPX vs. interp.')
        if len(_energy_list_dropped) <= 7:
            ax.legend()
        ax.grid()

        return ax

    def plot_tof_shape_interp(self, e_ev, source_to_detector_m, conv_proton, proton_params={},
                              t_interp=None, for_sum=False, logy=False, norm=False, marker='o', ax_mpl=None,
                              plt_arrow=True):
        """
        Plot each eV beam shape obtained from the fitting approach

        :param e_ev:
        :type e_ev:
        :param source_to_detector_m:
        :type source_to_detector_m:
        :param conv_proton:
        :type conv_proton:
        :param proton_params:
        :type proton_params:
        :param t_interp:
        :type t_interp:
        :param for_sum:
        :type for_sum:
        :param logy:
        :type logy:
        :param norm:
        :type norm:
        :param marker:
        :type marker:
        :param ax_mpl:
        :type ax_mpl:
        :param plt_arrow:
        :type plt_arrow:
        :return:
        :rtype:
        """
        if isinstance(e_ev, int) or isinstance(e_ev, float):
            e_ev = [e_ev]
        e_ev.sort()
        if t_interp is None:
            t_interp = self.t_us_mcnp
        self._make_shape(e_ev=e_ev, t_interp=t_interp, for_sum=for_sum, norm=norm,
                         source_to_detector_m=source_to_detector_m, print_tof=False,
                         conv_proton=conv_proton, proton_params=proton_params)
        _shape_tof_df_interp = self.shape_tof_df_interp
        _y_label = 'Flux (n/sterad/pulse)'

        if norm:
            _y_label = 'Ratio out of max flux of each energy'
        _x_tag = 'tof_us'
        _for_sum_s = ''
        if for_sum:
            _for_sum_s = ' for_sum'
        _norm_s = ''
        if norm:
            _norm_s = ' norm'
        _conv_proton_s = ''
        if conv_proton:
            _conv_proton_s = ' proton'
        _proton_param_s = ''
        for _param in proton_params.keys():
            _proton_param_s = _proton_param_s + '_' + _param + '_' + str(proton_params[_param])
        _details_s = _norm_s + _conv_proton_s + _proton_param_s + _for_sum_s
        _title_s = 'Energy dependent neutron pulse shape (interp.{})'.format(_details_s)

        if ax_mpl is None:
            fig, ax_mpl = plt.subplots()
        ax_mpl.set_prop_cycle(default_cycler)
        for each_e in e_ev:
            if not for_sum:
                _x_tag = str(each_e) + '_tof_us'
            if logy:
                ax_mpl.semilogy(_shape_tof_df_interp[_x_tag],
                                _shape_tof_df_interp[str(each_e)],
                                linestyle='--',
                                marker=marker,
                                mfc='none',
                                label=str(each_e) + ' eV (interp.{})'.format(_proton_param_s))
            else:
                ax_mpl.plot(_shape_tof_df_interp[_x_tag],
                            _shape_tof_df_interp[str(each_e)],
                            linestyle='--',
                            marker=marker,
                            mfc='none',
                            label=str(each_e) + ' eV (interp.{})'.format(_proton_param_s))
        if len(e_ev) <= 7:
            ax_mpl.legend()
        else:
            if plt_arrow:
                _plot_ev_arrow_as_legend(ax=ax_mpl, ev_list=e_ev)
        ax_mpl.set_ylabel(_y_label)
        ax_mpl.set_xlabel(u'Time (\u03BCs)')
        ax_mpl.grid()
        ax_mpl.set_title(_title_s)
        return ax_mpl

    def plot_proton_conv(self, e_ev, source_to_detector_m, conv_proton, sigma_list, tof=True,
                         t_interp=None, for_sum=False, logy=False, norm=False, plt_arrow=True):
        _list_of_dicts = _break_proton_param_list_to_dict(proton_param_list=sigma_list, proton_param_name='sigma')
        fig, ax1 = plt.subplots()
        for i, _e_param_dict in enumerate(_list_of_dicts):
            if tof:
                ax1 = self.plot_tof_shape_interp(e_ev=e_ev, source_to_detector_m=source_to_detector_m,
                                                 conv_proton=conv_proton, proton_params=_e_param_dict,
                                                 t_interp=t_interp, for_sum=for_sum, logy=logy, norm=norm, ax_mpl=ax1,
                                                 plt_arrow=plt_arrow, marker=marker_list[i])
            else:
                ax1 = self.plot_shape_interp(e_ev=e_ev, source_to_detector_m=source_to_detector_m,
                                             conv_proton=conv_proton, proton_params=_e_param_dict,
                                             t_interp=t_interp, for_sum=for_sum, logy=logy, norm=norm, ax_mpl=ax1,
                                             plt_arrow=plt_arrow, marker=marker_list[i])
        _for_sum_s = ''
        if for_sum:
            _for_sum_s = ' for_sum'
        _title_s = 'Energy dependent neutron pulse shape (interp.{}) with proton convolved'.format(_for_sum_s)
        ax1.set_title(_title_s)

        return ax1

    def make_shape(self, e_ev, source_to_detector_m, conv_proton, proton_params={},
                   t_interp=None, for_sum=False, norm=False, overwrite_csv=False):
        assert self.linear_df is not None
        assert self.model is not None
        if isinstance(e_ev, int) or isinstance(e_ev, float):
            e_ev = [e_ev]
        e_ev.sort()
        if t_interp is None:
            t_interp = self.t_us_mcnp
            # t_interp = self.t_us
        if isinstance(t_interp, int) or isinstance(t_interp, float):
            raise ValueError("'t_interp' must be a list or array.")
        t_interp.sort()

        # construct the str for the file name
        _distance_s = '_' + str(source_to_detector_m) + 'm'
        _for_sum_s = ''
        if for_sum:
            _for_sum_s = '_for_sum'
        _norm_s = ''
        if norm:
            _norm_s = '_norm'
        _conv_proton_s = ''
        if conv_proton:
            _conv_proton_s = '_proton'
        _proton_param_s = ''
        for _param in proton_params.keys():
            _proton_param_s = _proton_param_s + '_' + _param + '_' + str(proton_params[_param])

        _e_min = e_ev[0]
        _e_max = e_ev[-1]
        _e_nbr = len(e_ev) - 1
        _e_step = (_e_max - _e_min) / _e_nbr

        _e_str = '_eV_' + str(_e_min) + '_' + str(_e_max) + '_' + str(_e_step)

        _t_min = t_min_us
        _t_max = t_max_us
        _t_step = t_step_us

        _t_str = '_us_' + str(_t_min) + '_' + str(_t_max) + '_' + str(_t_step)

        assert self.model_used is not None
        _model_s = '_' + self.model_used + '.csv'

        _filename = 'TOF_shape' + _e_str + _t_str + _norm_s + _conv_proton_s + _proton_param_s \
                    + _for_sum_s + _distance_s + _model_s
        _shape_tof_df_dir = os.path.join(self.result_neutron_folder, _filename)
        self.shape_tof_df_dir = _shape_tof_df_dir

        # File exists
        if os.path.isfile(_shape_tof_df_dir):
            print("\u2705 '{}' exists...".format(_shape_tof_df_dir))
            if overwrite_csv:
                # Override==True, perform making shape and overwrite the .csv file
                print("File overwriting...")
                print("New beam shape generation starts...")
                # Making starts
                self._make_shape(e_ev=e_ev, t_interp=t_interp, for_sum=for_sum, norm=norm,
                                 source_to_detector_m=source_to_detector_m,
                                 save_dir=_shape_tof_df_dir,
                                 print_tof=True,
                                 conv_proton=conv_proton,
                                 proton_params=proton_params,
                                 )
                print("File overwritten.")
            else:
                # Override==False, read the .csv file
                self.shape_tof_df_interp = pd.read_csv(_shape_tof_df_dir)
                print("TOF neutron beam shape file loaded.")
                self.proton_pulse.make_new_shape(proton_params=proton_params)  # Making sure proton shape params updated

        # File not exists, perform fitting
        else:
            print("\u274C No previous TOF neutron beam shape file named '{}' detected.".format(_shape_tof_df_dir))
            print("Beam shape generation starts...")
            # Making starts
            self._make_shape(e_ev=e_ev, t_interp=t_interp, for_sum=for_sum, norm=norm,
                             source_to_detector_m=source_to_detector_m,
                             save_dir=_shape_tof_df_dir,
                             print_tof=True,
                             conv_proton=conv_proton,
                             proton_params=proton_params,
                             )

    def _make_shape(self, e_ev, t_interp, for_sum, norm, source_to_detector_m, print_tof,
                    conv_proton, proton_params={}, save_dir=None):
        assert self.linear_df is not None
        assert self.model is not None
        if isinstance(e_ev, int) or isinstance(e_ev, float):
            e_ev = [e_ev]
        e_ev.sort()
        if t_interp is not None:
            t_interp.sort()
        else:
            t_interp = self.t_us_mcnp
        if isinstance(t_interp, int) or isinstance(t_interp, float):
            raise ValueError("'t_interp' must be a list or array.")
        # _t_shift = 0
        # if conv_proton:
        #     _t_shift = self.proton_pulse._t_shift
        #     # t used to shift convolution from 'full' to 'same'

        _param_df_interp = self._interpolate_param(e_ev=e_ev).set_index('E_eV')
        _shape_df_interp = pd.DataFrame()
        _shape_df_interp['t_us'] = t_interp
        # _shape_df_interp['t_us'] = t_interp - _t_shift
        _shape_tof_df_interp = pd.DataFrame()
        _tof_us_dict = {}
        _tof_total_us_array = []

        self.proton_pulse.make_new_shape(proton_params=proton_params)  # Making sure proton shape params updated

        if print_tof:
            print('For {} (m)'.format(source_to_detector_m))

        for _each_e in e_ev:
            _t_us, _array = self._make_single_shape(e_ev=_each_e,
                                                    t_us=t_interp,
                                                    param_df=_param_df_interp,
                                                    conv_proton=conv_proton,
                                                    )
            if not norm:
                _array = _array * _param_df_interp['f_max'][_each_e]
            _array[_array < 0] = 0
            _array = _array.round(5)
            _shape_df_interp[_each_e] = _array
            _tof_diff_us = ev_to_s(offset_us=0, source_to_detector_m=source_to_detector_m, array=_each_e) * 1e6
            if print_tof:
                print('{} (eV) neutron spent {} (us)'.format(_each_e, _tof_diff_us))
            _tof_us_dict[_each_e] = _tof_diff_us
            _current_tof_us = t_interp + _tof_diff_us

            _tof_total_us_array = np.append(_tof_total_us_array, _current_tof_us)

            if not for_sum:
                # _shape_tof_df_interp[str(_each_e) + '_tof_us'] = _current_tof_us - _t_shift
                _shape_tof_df_interp[str(_each_e) + '_tof_us'] = _t_us + _tof_diff_us
                _shape_tof_df_interp[str(_each_e)] = _array

        self.shape_df_interp = _shape_df_interp
        self.tof_us_dict = _tof_us_dict

        _tof_total_us_array.sort()  # list of all time that exist in all energy

        if for_sum:
            # _shape_tof_df_interp['tof_us'] = _tof_total_us_array - _t_shift
            if print_tof is True:
                print('Making shape for:')
            for _each_e in e_ev:
                __tof_diff_us = _tof_us_dict[_each_e]
                _current_t_without_tof = _tof_total_us_array - __tof_diff_us
                if print_tof is True:
                    print('{} (eV) neutron ...'.format(_each_e))
                _t_us, _array = self._make_single_shape(e_ev=_each_e,
                                                        t_us=_current_t_without_tof,
                                                        param_df=_param_df_interp,
                                                        conv_proton=conv_proton,
                                                        )
                if not norm:
                    _array = _array * _param_df_interp['f_max'][_each_e]
                _array[_array < 0] = 0
                _array = _array.round(5)
                _shape_tof_df_interp['tof_us'] = _t_us + __tof_diff_us
                _shape_tof_df_interp[str(_each_e)] = _array

        self.shape_tof_df_interp = _shape_tof_df_interp

        # Save shape_tof_df_interp as .csv
        if save_dir is not None:
            self.shape_tof_df_interp.to_csv(save_dir, index=False)
            print("TOF neutron beam shape file has been saved at '{}'".format(save_dir))

    def _make_single_shape(self, e_ev, t_us, param_df, conv_proton):
        # if not isinstance(e_ev, int) or isinstance(e_ev, float):
        #     raise ValueError("'e_ev' must be a number for single shape generation.")
        if isinstance(t_us, int) or isinstance(t_us, float):
            raise ValueError("'t_us' must be a list or array for shape interpolation.")
        # if e_ev < 1 or e_ev > 500:
        #     raise ValueError()
        t_us = np.array(t_us)
        _my_model = self.model
        for _each_param in self.model_param_names:
            _my_model.set_param_hint(_each_param, value=param_df[_each_param][e_ev])
        _params = _my_model.make_params()

        _t_shift = 0
        if not conv_proton:
            _array = _my_model.eval(_params, t=t_us)  # lmfit.model.eval() returns np.ndarray
        else:
            _t_shift = self.proton_pulse._t_shift
            # t used to shift convolution from 'full' to 'same'
            _array_for_conv_proton = _my_model.eval(_params, t=self.t_us_conv_proton)
            _proton_x = np.array(self.proton_pulse.new_shape_df['t_ns'] / 1e3 + self.t_us_conv_proton[-1])
            _proton_y = np.array(self.proton_pulse.new_shape_df['intensity'])
            _conv_y = np.convolve(_array_for_conv_proton, _proton_y, mode='full')
            _conv_x = np.append(self.t_us_conv_proton, _proton_x[1:])
            _array_function = interp1d(x=_conv_x, y=_conv_y, kind='cubic', bounds_error=False, fill_value=0)
            _array = _array_function(t_us)
            # print(t_us)
            # print(len(t_us))

        assert len(t_us) == len(_array)
        return t_us - _t_shift, _array

    def _interpolate_param(self, e_ev):

        _linear_df = self.linear_df.set_index('param_name')
        _param_df_interp = pd.DataFrame()
        _param_df_interp['E_eV'] = e_ev

        for each_param in _linear_df.index:
            _param_df_interp[each_param] = loglog_linear(x=e_ev,
                                                         slope=_linear_df['slope'][each_param],
                                                         intercept=_linear_df['intercept'][each_param])
        return _param_df_interp

    def plot_params_vs_e(self, loglog=True):
        assert self.param_df is not None
        self.param_df.set_index('E_eV').plot(loglog=loglog, style='.')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Fitted parameter value')

    def fit_params(self, show_init=True, check_each=False, save_fig=False, overwrite_csv=False, loglog_fit=True):
        if self.param_df is None:
            raise ValueError("'NeutronPulse._fit_shape_proton()' must be applied before 'NeutronPulse.fit_params'")

        _e_min = str(self.e_min) + 'eV_'
        _e_max = str(self.e_max) + 'eV_'
        _model_s = self.model_used + '.csv'
        _filename = 'Loglog_linear_within_' + _e_min + _e_max + _model_s
        self.linear_df_dir = os.path.join(self.result_neutron_folder, _filename)

        # File exists
        if os.path.isfile(self.linear_df_dir):
            print("\u2705 '{}' exists...".format(self.linear_df_dir))
            if overwrite_csv:
                # Override==True, perform fitting and overwrite the .csv file
                print("File overwriting...")
                print("New fitting starts...")
                # Fitting starts
                self._fit_params(show_init=show_init,
                                 check_each=check_each,
                                 save_fig=save_fig,
                                 loglog_fit=loglog_fit)
                print("File overwritten.")
            else:
                # Override==False, read the .csv file
                self.linear_df = pd.read_csv(self.linear_df_dir)
                print("Parameters linear fitted file loaded.")

        # File not exists, perform fitting
        else:
            print("\u274C No previous fitting file detected.\nNew fitting starts...")
            # Fitting starts
            self._fit_params(show_init=show_init,
                             check_each=check_each,
                             save_fig=save_fig,
                             loglog_fit=loglog_fit)

    def _fit_params(self, show_init, check_each, save_fig, loglog_fit):

        _param_df = self.param_df
        e_log = np.log10(_param_df['E_eV'])
        param_name_list = list(_param_df.columns.drop('E_eV'))
        linear_dict = {}
        for each_param in param_name_list:
            linear_dict[each_param] = self.__fit_params(y=_param_df[each_param],
                                                        x=_param_df['E_eV'],
                                                        x_log=e_log,
                                                        name=each_param,
                                                        check_each=check_each,
                                                        show_init=show_init,
                                                        save_fig=save_fig,
                                                        loglog_fit=loglog_fit
                                                        )
        # self.linear_dict = linear_dict
        self.linear_df = self._form_linear_df(linear_dict=linear_dict, save=True)

    def _form_linear_df(self, linear_dict, save=False):
        _key_list = list(linear_dict.keys())
        linear_param_list = list(linear_dict[_key_list[0]].keys())

        _dict = {}
        for each_param in linear_param_list:
            _dict[each_param] = []
        for _key in _key_list:
            for each_param in linear_param_list:
                _dict[each_param].append(linear_dict[_key][each_param])
        _df = pd.DataFrame()
        for each_param in linear_param_list:
            _df[each_param] = _dict[each_param]
        _df.insert(0, 'param_name', _key_list)

        if save is True:
            assert self.result_neutron_folder is not None
            assert self.linear_df_dir is not None
            _df.to_csv(path_or_buf=self.linear_df_dir, index=False)
            print("Parameters linear fitting file has been saved at '{}'".format(self.linear_df_dir))

        return _df

    def __fit_params(self, y, x, x_log, name, show_init, check_each, save_fig, loglog_fit):

        # Guess to make param
        y_log = np.log10(y)
        _temp_guess_model = LinearModel()
        params = _temp_guess_model.guess(y_log, x_log)

        # Load actual loglog linear model
        my_param_model = Model(loglog_linear)
        assert my_param_model.independent_vars == ['x']
        assert my_param_model.param_names == ['slope', 'intercept']

        # Fit the model (NOTE: lmfit leastsq gives different results!!!???)
        if loglog_fit:
            # fit in loglog space
            out = _temp_guess_model.fit(y_log, params, x=x_log)
            _x_label = 'Log E (eV)'
            _y_label = 'Log value (arb. unit)'
        else:
            # fit in real space
            out = my_param_model.fit(y, params, x=x)
            _y_label = 'Parameter value (arb. unit)'
            _x_label = 'E (eV)'

        if check_each:
            # assert self.result_shape_fit is not None
            if name == 'f_max':
                if loglog_fit:
                    _y_label = 'Log flux (neutrons/s/cm2)'
                else:
                    _y_label = 'Flux (neutrons/s/cm2)'

            if show_init is True:
                out.plot(xlabel=_x_label, ylabel=_y_label)
            else:
                out.plot(xlabel=_x_label, ylabel=_y_label, initfmt='None')

            plt.title(name, fontsize=12)
            if loglog_fit is False:
                plt.yscale('log')
                plt.xscale('log')
                plt.xlim(left=self.e_min * 0.8, right=self.e_max * 1.2)
                plt.ylim(bottom=np.amin(y) * 0.8, top=np.amax(y) * 1.2)

            if save_fig:
                # Check and make dir to save
                assert self.result_neutron_folder is not None
                assert self.model_used is not None
                _model_s = self.model_used + '_' + name + '.png'
                _filename = 'Neutron_pulse_param_fit_' + _model_s
                _dir_to_save = os.path.join(self.result_neutron_folder, _filename)
                plt.savefig(_dir_to_save, dpi=300, transparent=False)
                plt.close()
            else:
                plt.show()
        elif save_fig:
            raise ValueError("'check_each' has to be 'True' in order to save figure")

        return out.best_values

    def __set_model(self, model_index):
        if model_index == 1:
            my_model = Model(ikeda_carpenter)
        elif model_index == 2:
            my_model = Model(cole_windsor)
        elif model_index == 3:
            my_model = Model(pseudo_voigt)
        elif model_index == 4:
            my_model = Model(ikeda_carpenter_jparc)
        elif model_index == 5:
            my_model = Model(cole_windsor_jparc)
        else:
            raise ValueError("Model index not exists, please refer to: '{}' ".format(self.model_map))

        self.model = my_model
        self.model_index = model_index
        self.model_used = self.model_map[model_index]
        self.model_param_names = my_model.param_names

    def fit_shape(self, e_min, e_max,
                  drop=False, norm=True,
                  show_init=True, check_each=False, save_fig=False, overwrite_csv=False):
        # [1: 'ikeda_carpenter', 2: 'cole_windsor', 3: 'pseudo_voigt']
        self.e_min = e_min
        self.e_max = e_max
        self._form_energy_list_dropped(e_min=e_min, e_max=e_max)
        # File name of fitted_param_df.csv
        assert self.model_used is not None
        _e_min = str(self.e_min) + 'eV_'
        _e_max = str(self.e_max) + 'eV_'
        _model_s = self.model_used + '.csv'
        _filename = 'Neutron_fitted_params_' + _e_min + _e_max + _model_s
        self.param_df_dir = os.path.join(self.result_neutron_folder, _filename)

        # File exists
        if os.path.isfile(self.param_df_dir):
            print("\u2705 '{}' exists...".format(self.param_df_dir))
            if overwrite_csv:
                # Override==True, perform fitting and overwrite the .csv file
                print("File overwriting...")
                print("New fitting starts...")
                # Fitting starts
                self._fit_shape_params(drop=drop, norm=norm, show_init=show_init, check_each=check_each,
                                       save_fig=save_fig)
                print("File overwritten.")
            else:
                # Override==False, read the .csv file
                self.param_df = pd.read_csv(self.param_df_dir)
                print("Fitted parameters file loaded.")

        # File not exists, perform fitting
        else:
            print("\u274C No previous fitting file detected.\nNew fitting starts...")
            # Fitting starts
            self._fit_shape_params(drop=drop, norm=norm, show_init=show_init, check_each=check_each, save_fig=save_fig)

    def _fit_shape_params(self, drop, norm, show_init, check_each, save_fig):
        # Fitting starts
        param_dict_fitted = {}
        for each_e in self.shape_dict_mcnp.keys():
            if self.e_min <= each_e <= self.e_max:
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
                f = self.shape_dict_mcnp[each_e][_data_used][_flux_used]
                t = self.shape_dict_mcnp[each_e][_data_used]['t_us']
                param_dict_fitted[each_e]['fitted_params'] = self.__fit_shape_params(f=f,
                                                                                     t=t,
                                                                                     e=each_e,
                                                                                     model_index=self.model_index,
                                                                                     check_each=check_each,
                                                                                     show_init=show_init,
                                                                                     save_fig=save_fig)
        # Organize fitted parameters into pd.DataFrame
        self.param_df = self._form_params_df(param_dict=param_dict_fitted, save=True)

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
            f_max_list.append(self.shape_dict_mcnp[each_e]['f_max'])
            for each_param in param_list_fitted:
                _dict[each_param].append(param_dict[each_e]['fitted_params'][each_param])
        _df = pd.DataFrame()
        _df['E_eV'] = e_list_fitted
        for each_param in param_list_fitted:
            _df[each_param] = _dict[each_param]
        _df['f_max'] = f_max_list

        # print("------------------------------------------------------------------------------------------")
        # print("NOTE: 'f_max' in the params_df is NOT a fitted parameter, it is the maximum in each shape.")
        # print("------------------------------------------------------------------------------------------")
        # print("                                    Fitted Parameter                                      ")
        # print("------------------------------------------------------------------------------------------")
        # print(_df)

        if save is True:
            assert self.result_neutron_folder is not None
            assert self.param_df_dir is not None
            _df.to_csv(path_or_buf=self.param_df_dir, index=False)
            print("Parameters fitting file ({}) has been saved at '{}'".format(self.model_used, self.param_df_dir))

        return _df

    def __fit_shape_params(self, f, t, e, model_index, show_init, check_each, save_fig):
        _model_map = self.model_map
        verbose = False

        if model_index != self.model_index:
            self.__set_model(model_index)
        my_model = self.model
        # ikeda_carpenter
        if model_index == 1:
            assert my_model.param_names == ['alpha', 'beta', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']

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
            assert my_model.param_names == ['alpha', 'beta', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']

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
            assert my_model.param_names == ['sig1', 'sig2', 'gamma', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']

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
            assert my_model.param_names == ['sig1', 'sig2', 'gam1', 'gam2', 'fraction', 't0', 'norm_factor']
            assert my_model.independent_vars == ['t']

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
            _x_label = 'Time (\u03BCs)'
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

    def _export_total(self, filename=None):
        assert self.shape_total_df is not None

        if filename is None:
            self.shape_total_df.to_clipboard(excel=True)
        else:
            self.shape_total_df.to_csv(filename, index=False)

    def _export_each_to_csv(self):
        for index, each_energy in enumerate(self.shape_dict_mcnp.keys()):
            df = self.shape_dict_mcnp[each_energy]
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
        self._shape_df = _load_proton_pulse(path)
        self._shape_df_fit = None
        self.params = None
        self.new_params = None
        self.new_shape_df = None
        self.model = None
        self._t_shift = None
        self.result = self._fit_shape_proton(print_params=True)

    def __get_t_shift(self):
        self._t_shift = round(self.new_shape_df['t_ns'].median() / 10) / 100

    def _fit_shape_proton(self, print_params=False):
        t_ns = self._shape_df['t_ns']
        intensity = self._shape_df['intensity']
        _model = lmfit.models.GaussianModel()
        _proton_params = _model.guess(data=intensity, x=t_ns)
        result = _model.fit(data=intensity, x=t_ns, params=_proton_params)
        self.model = _model
        self.params = result.params
        sig_fit_s = str(round(self.params.valuesdict()['sigma'], 2))
        # if self._shape_df_fit is None:
        self._shape_df_fit = pd.DataFrame()
        self._shape_df_fit['t_ns'] = self._shape_df['t_ns']
        _new_array = self.model.eval(params=result.params, x=self._shape_df['t_ns'])
        self._shape_df_fit['Fit (\u03C3={})'.format(sig_fit_s)] = _new_array.round(5)

        _new_shape_df = pd.DataFrame()
        _new_shape_df['t_ns'] = self._shape_df_fit['t_ns']
        _new_shape_df['intensity'] = _new_array.round(5)
        self.new_shape_df = _new_shape_df
        self.__get_t_shift()

        if print_params:
            print("Fitted params for raw proton pulse shape:")
            result.params.pretty_print()
        return result

    def make_new_shape(self, proton_params={}, print_params=False):
        _check_proton_params_dict_keys(proton_params)
        assert self.model is not None
        _params_guess = self.model.guess(data=self._shape_df['intensity'], x=self._shape_df['t_ns'])
        result = self.model.fit(data=self._shape_df['intensity'], x=self._shape_df['t_ns'], params=_params_guess)
        _new_params = result.params
        for each_par in proton_params.keys():
            _new_params.add(each_par, proton_params[each_par])
        if print_params:
            print("Fitted params for simulated proton pulse shape:")
            _new_params.pretty_print()
        sig_fit_s = str(round(_new_params.valuesdict()['sigma'], 2))
        _new_shape_df = pd.DataFrame()
        _new_shape_df['t_ns'] = self._shape_df['t_ns']
        _new_array = self.model.eval(params=_new_params, x=self._shape_df['t_ns'])
        _new_shape_df['intensity'] = _new_array.round(5)
        self._shape_df_fit['Fit (\u03C3={})'.format(sig_fit_s)] = _new_array.round(5)
        self.new_params = _new_params
        self.new_shape_df = _new_shape_df
        self.__get_t_shift()
        return _new_shape_df

    def plot(self):
        fig, ax = plt.subplots()
        # sig1 = str(round(self.params.valuesdict()['sigma'], 2))
        ax.plot(self._shape_df['t_ns'], self._shape_df['intensity'], label='Raw')
        # ax.plot(_x, self._shape_df_fit['intensity'], 'b:', label='Fit (\u03C3={})'.format(sig1))
        if self._shape_df_fit is not None:
            # sig2 = str(round(self.new_params.valuesdict()['sigma'], 2))
            # ax.plot(_x, self.new_shape_df['intensity'], 'r--', label='Fit (\u03C3={})'.format(sig2))
            ax = self._shape_df_fit.set_index('t_ns').plot(ax=ax)
        ax.set_title('Proton pulse raw vs. fitting')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Intensity (a.u.)')
        ax.legend()
        ax.grid()
        return ax


def _check_proton_params_dict_keys(proton_params):
    _key_list = ['sigma', 'center', 'amplitude']
    for each_par in proton_params.keys():
        if each_par not in _key_list:
            raise ValueError("'{}' is not a valid param name, accepted keys are: '{}'.".format(each_par, _key_list))


def _break_proton_param_list_to_dict(proton_param_list, proton_param_name):
    _list_of_dicts = []
    for _each in proton_param_list:
        _list_of_dicts.append({proton_param_name: _each})
    return _list_of_dicts


# Functions to load files #
def _load_neutron_total_shape(path):
    p = np.genfromtxt(path)
    pp = p.T
    df1 = pd.DataFrame()
    for i, col in enumerate(pp):
        df1[i] = col
    col_name_1 = ['E_eV', 'l_angstrom', 'b(E)', 'bs(E)', 'b(l)', 'bs(l)']
    df1.columns = col_name_1
    df1['f(E)'] = df1['b(E)'] * df1['E_eV']
    df1['s(E)'] = df1['bs(E)'] * df1['E_eV']
    df1['f(l)'] = df1['b(l)'] * df1['l_angstrom']
    df1['s(l)'] = df1['bs(l)'] * df1['l_angstrom']

    return df1


def _load_neutron_each_shape(path, export=False):
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
        b = []
        bs = []
        f = []
        s = []
        df = pd.DataFrame()
        for each_line in q:
            if each_energy == each_line[1]:
                t_us.append(each_line[0])
                e_ev.append(each_line[1])
                b.append(each_line[2])
                f.append(each_line[2] * each_line[1])  # brightness * E
                bs.append(each_line[3])
                s.append(each_line[3] * each_line[1])  # brightness s * E
        f_max = np.amax(f)
        df['t_us'] = t_us
        df['b'] = b
        df['bs'] = bs
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


def _shape_dict_to_dfs(shape_dict, t_max=None):
    _first_key = list(shape_dict.keys())[0]
    _t_us = shape_dict[_first_key]['data_raw']['t_us']
    _shape_df = pd.DataFrame()
    _shape_df_norm = pd.DataFrame()
    _shape_df['t_us'] = _t_us
    _shape_df_norm['t_us'] = _t_us
    for each_key in shape_dict.keys():
        _shape_df[each_key] = shape_dict[each_key]['data_raw']['f']
        _shape_df_norm[each_key] = shape_dict[each_key]['data_raw']['f_norm']
    # _shape_df = _shape_df.drop(_shape_df[_shape_df.t_us < t_min_us].index)
    # _shape_df.reset_index(drop=True, inplace=True)
    # _shape_df = _shape_df.drop(_shape_df[_shape_df.t_us < t_min_us].index)
    # _shape_df.reset_index(drop=True, inplace=True)
    # , t_min_us, t_max_us
    if t_max is not None:
        _shape_df = _shape_df[_shape_df.t_us <= t_max]
        _shape_df_norm = _shape_df_norm[_shape_df_norm.t_us <= t_max]
    return _shape_df, _shape_df_norm


def _load_proton_pulse(path=proton_path):
    _df = pd.read_csv(path, sep=' ', skiprows=1, header=None)
    _df.columns = ['t_ns', 'intensity']
    df = _df.round(5)
    df[df < 0] = 0
    return df


def _plot_ev_arrow_as_legend(ax, ev_list):
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0.5, alpha=0.4)
    ax.annotate(str(ev_list[-1]) + ' eV',
                xy=(.65, .65), xycoords='figure fraction',
                xytext=(.21, .83), textcoords='figure fraction',
                arrowprops=dict(arrowstyle="fancy",
                                # linestyle="dashed",
                                color="k",
                                patchB=None,
                                shrinkB=0,
                                connectionstyle="arc3",
                                alpha=0.55
                                ),
                bbox=bbox_props
                )
    ax.annotate(str(ev_list[0]) + ' eV',
                xy=(.63, .62), xycoords='figure fraction',
                bbox=bbox_props
                )
    return ax
