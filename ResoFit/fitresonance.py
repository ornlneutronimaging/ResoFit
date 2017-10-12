import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation
import numpy as np
from lmfit import minimize
import re
from ResoFit._gap_functions import y_gap_for_fitting
from ResoFit._gap_functions import y_gap_for_iso_fitting
import periodictable as pt
from ResoFit._utilities import Layer
import pprint


class FitResonance(Experiment):
    fit_result = None
    fitted_density_gcm3 = None
    fitted_thickness_mm = None
    fitted_residual = None
    fitted_gap = None
    fitted_fjac = None
    fitted_layer = None
    fitted_simulation = None
    layer_list = None
    raw_layer = None
    fitted_iso_result = None
    fitted_iso_residual = None
    isotope_stack = {}

    def __init__(self, spectra_file, data_file,
                 calibrated_offset_us, calibrated_source_to_detector_m,
                 folder, repeat=1, baseline=False,
                 norm_to_file=None, slice_start=None, slice_end=None,
                 energy_min=1e-5, energy_max=1000, energy_step=0.01):
        super().__init__(spectra_file=spectra_file, data_file=data_file, folder=folder, repeat=repeat)
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.energy_step = energy_step
        self.calibrated_offset_us = calibrated_offset_us
        self.calibrated_source_to_detector_m = calibrated_source_to_detector_m
        self.raw_layer = None
        self.slice(slice_start=slice_start, slice_end=slice_end)
        self.baseline = baseline
        if norm_to_file is not None:
            self.norm_to(norm_to_file)
        self.exp_x_interp, self.exp_y_interp = self.xy_scaled(energy_min=self.energy_min,
                                                              energy_max=self.energy_max,
                                                              energy_step=self.energy_step,
                                                              angstrom=False, transmission=False,
                                                              offset_us=self.calibrated_offset_us,
                                                              source_to_detector_m=self.calibrated_source_to_detector_m,
                                                              baseline=self.baseline)

    def fit(self, raw_layer, vary='density', isotope=False, isotope_vary_tag_list=[], each_step=False):
        if vary not in ['density', 'thickness', 'none']:
            raise ValueError("'vary=' can only be one of ['density', 'thickness', 'none']")
        # Default vary is: 'density'
        thickness_vary_tag = False
        density_vary_tag = True
        if vary == 'thickness':
            thickness_vary_tag = True
            density_vary_tag = False
        if vary == 'none':
            density_vary_tag = False
        self.raw_layer = raw_layer
        '''Load params'''
        self.layer_list = list(raw_layer.info.keys())
        params_for_fit = Parameters()
        for _each_layer in self.layer_list:
            if self.raw_layer.info[_each_layer]['density']['value'] is np.NaN:
                self.raw_layer.info[_each_layer]['density']['value'] = pt.elements.isotope(_each_layer).density
            params_for_fit.add('thickness_mm_' + _each_layer,
                               value=self.raw_layer.info[_each_layer]['thickness']['value'],
                               vary=thickness_vary_tag,
                               min=0)
            params_for_fit.add('density_gcm3_' + _each_layer,
                               value=self.raw_layer.info[_each_layer]['density']['value'],
                               vary=density_vary_tag,
                               min=0)
        # Print before
        print("Params before '{}' fitting:".format(vary))
        params_for_fit.pretty_print()
        # Fitting
        self.fit_result = minimize(y_gap_for_fitting, params_for_fit, method='leastsq',
                                   args=(self.exp_x_interp, self.exp_y_interp, self.layer_list,
                                         self.energy_min, self.energy_max, self.energy_step, each_step))
        # Print after
        print("Params after '{}' fitting:".format(vary))
        self.fit_result.__dict__['params'].pretty_print()
        # Print chi^2
        self.fitted_residual = self.fit_result.__dict__['residual']
        print("Fitting chi^2 : {}\n".format(sum(self.fitted_residual ** 2)))

        '''Export fitted params as Layer()'''
        # Save the fitted 'density' or 'thickness' in Layer()
        self.fitted_layer = Layer()
        for _each_layer in self.layer_list:
            self.fitted_layer.add_layer(layer=_each_layer,
                                        thickness_mm=self.fit_result.__dict__['params'].valuesdict()[
                                            'thickness_mm_' + _each_layer],
                                        density_gcm3=self.fit_result.__dict__['params'].valuesdict()[
                                            'density_gcm3_' + _each_layer])
        # self.fitted_fjac = self.fit_result.__dict__['fjac']
        # print(self.fit_result.__dict__['fjac'][0])

        '''Create fitted simulation'''
        self.fitted_simulation = Simulation(energy_min=self.energy_min,
                                            energy_max=self.energy_max,
                                            energy_step=self.energy_step)
        for each_layer in self.layer_list:
            self.fitted_simulation.add_layer(layer=each_layer,
                                             layer_thickness_mm=self.fitted_layer.info[each_layer]['thickness'][
                                                 'value'],
                                             layer_density_gcm3=self.fitted_layer.info[each_layer]['density']['value'])

        return self.fit_result

    def fit_iso(self, layer, each_step=False):
        params_for_iso_fit = Parameters()
        self.isotope_stack[layer] = {'list': self.fitted_simulation.o_reso.stack[layer][layer]['isotopes']['list'],
                                     'ratios': self.fitted_simulation.o_reso.stack[layer][layer]['isotopes'][
                                         'isotopic_ratio']}
        _formatted_isotope_list = []
        _params_name_list = []
        _restriction_list = []
        for _isotope_index in range(len(self.isotope_stack[layer]['list'])):
            _formatted_isotope_name = self.isotope_stack[layer]['list'][_isotope_index].replace('-', '_')
            _formatted_isotope_list.append(_formatted_isotope_name)
            _params_name_list.append('isotope_ratio_' + _formatted_isotope_name)

        for _name_index in range(len(_params_name_list)):
            params_for_iso_fit.add(_params_name_list[_name_index],
                                   value=self.isotope_stack[layer]['ratios'][_name_index],
                                   min=0,
                                   max=1)

        # Restriction is not working as expected ###
        # Create expr to restrict the total to be 1.0
        # for _each_param_name in _params_name_list:
        #     _params_name_list_temp = _params_name_list[:]
        #     _params_name_list_temp.remove(_each_param_name)
        #     _params_name_list_temp.insert(0, '1.0')
        #     _restriction_list.append(' - '.join(_params_name_list_temp))
        #     print(_params_name_list_temp)
        # print(_restriction_list)
        # params_for_iso_fit.pretty_print()
        # for _i in range(len(_restriction_list)):
        #     params_for_iso_fit[_params_name_list[_i]].set(expr=_restriction_list[_i], vary=True)

        # Print params before
        print("Params before 'isotope' fitting:")
        params_for_iso_fit.pretty_print()
        # Fitting
        self.fitted_iso_result = minimize(y_gap_for_iso_fitting, params_for_iso_fit, method='leastsq',
                                          args=(self.exp_x_interp, self.exp_y_interp, layer, _formatted_isotope_list,
                                                self.fitted_simulation, each_step))
        # Print params after
        print("Params after 'isotope' fitting:")
        self.fitted_iso_result.__dict__['params'].pretty_print()
        # Print chi^2
        self.fitted_iso_residual = self.fitted_iso_result.__dict__['residual']
        print("Fit iso chi^2 : {}\n".format(sum(self.fitted_iso_residual ** 2)))

        return

    def molar_conc(self):
        molar_conc_units = 'mol/cm3'
        print("Molar-conc. ({})\tBefore_fit\tAfter_fit".format(molar_conc_units))
        for _each_layer in self.layer_list:
            molar_mass_value = self.fitted_simulation.o_reso.stack[_each_layer][_each_layer]['molar_mass']['value']
            molar_mass_units = self.fitted_simulation.o_reso.stack[_each_layer][_each_layer]['molar_mass']['units']
            # Adding molar_mass to fitted_layer info
            self.fitted_layer.info[_each_layer]['molar_mass']['value'] = molar_mass_value
            self.fitted_layer.info[_each_layer]['molar_mass']['units'] = molar_mass_units
            # Adding molar_mass to raw_layer info
            self.raw_layer.info[_each_layer]['molar_mass']['value'] = molar_mass_value
            self.raw_layer.info[_each_layer]['molar_mass']['units'] = molar_mass_units
            # Adding molar_concentration to fitted_layer info
            molar_conc_value = self.fitted_layer.info[_each_layer]['density']['value'] / molar_mass_value
            self.fitted_layer.info[_each_layer]['molar_conc']['value'] = molar_conc_value
            self.fitted_layer.info[_each_layer]['molar_conc']['units'] = molar_conc_units
            # Calculate starting molar_concentration and fitted_layer info
            start_molar_conc_value = self.raw_layer.info[_each_layer]['density']['value'] / molar_mass_value
            self.raw_layer.info[_each_layer]['molar_conc']['value'] = start_molar_conc_value
            self.raw_layer.info[_each_layer]['molar_conc']['units'] = molar_conc_units
            # molar_conc_output[_each_layer] = {'Before_fit': start_molar_conc_value,
            #                                   'After_fit': molar_conc_value}
            print("{}\t{}\t{}".format(_each_layer, start_molar_conc_value, molar_conc_value))

        return self.fitted_layer.info

    def plot_before(self):
        # Form signals from raw raw_layer
        simulation = Simulation(energy_min=self.energy_min,
                                energy_max=self.energy_max,
                                energy_step=self.energy_step)
        for each_layer in self.layer_list:
            simulation.add_layer(layer=each_layer,
                                 layer_thickness_mm=self.raw_layer.info[each_layer]['thickness']['value'],
                                 layer_density_gcm3=self.raw_layer.info[each_layer]['density']['value'])
        simu_x, simu_y = simulation.xy_simu(angstrom=False, transmission=False)

        # Get plot labels
        simu_label = 'Ideal'
        exp_label = 'Exp'
        exp_interp_label = 'exp_interp'
        for each_layer in self.layer_list:
            simu_label = simu_label + '_' + each_layer
            exp_label = exp_label + '_' + each_layer
            exp_interp_label = exp_interp_label + '_' + each_layer

        # Plot
        plt.plot(simu_x, simu_y,
                 'b-', label=simu_label, markersize=1)
        plt.plot(self.x_raw(angstrom=False, offset_us=self.calibrated_offset_us,
                            source_to_detector_m=self.source_to_detector_m),
                 self.y_raw(transmission=False, baseline=self.baseline),
                 'ro', label=exp_label, markersize=1)

        plt.title('Before fitting')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Attenuation')
        plt.ylim(ymax=1.01)
        plt.xlim(0, self.energy_max)
        plt.legend(loc='best')
        plt.show()

    def plot_after(self, error=True):
        # Form signals from fitted raw_layer
        if self.fitted_simulation is None:
            self.fitted_simulation = Simulation(energy_min=self.energy_min,
                                                energy_max=self.energy_max,
                                                energy_step=self.energy_step)
            for each_layer in self.layer_list:
                self.fitted_simulation.add_layer(layer=each_layer,
                                                 layer_thickness_mm=self.fitted_layer.info[each_layer]['thickness'][
                                                     'value'],
                                                 layer_density_gcm3=self.fitted_layer.info[each_layer]['density'][
                                                     'value'])
        simu_x, simu_y = self.fitted_simulation.xy_simu(angstrom=False, transmission=False)

        # Get plot labels
        simu_label = 'Ideal'
        exp_label = 'Exp'
        exp_interp_label = 'Exp_interp'
        for each_layer in self.layer_list:
            simu_label = simu_label + '_' + each_layer
            exp_label = exp_label + '_' + each_layer
            exp_interp_label = exp_interp_label + '_' + each_layer
        plt.plot(simu_x, simu_y,
                 'b-', label=simu_label, markersize=1)

        plt.plot(self.x_raw(angstrom=False, offset_us=self.calibrated_offset_us,
                            source_to_detector_m=self.source_to_detector_m),
                 self.y_raw(transmission=False, baseline=self.baseline),
                 'ro', label=exp_label, markersize=1)

        if error is True:
            # Plot fitting differences
            plt.plot(simu_x, self.fitted_residual - 0.2, 'g-', label='Diff.', alpha=0.8)

        plt.title('Best fit')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Attenuation')
        plt.ylim(ymax=1.01)
        plt.xlim(0, self.energy_max)
        plt.legend(loc='best')
        plt.show()
