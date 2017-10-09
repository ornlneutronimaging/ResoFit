import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation
import numpy as np
from lmfit import minimize
import re
from ResoFit._gap_functions import y_gap_for_fitting
import periodictable as pt
from ResoFit._utilities import Layer


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

    def __init__(self, spectra_file, data_file,
                 calibrated_offset_us, calibrated_source_to_detector_m,
                 folder='data', repeat=1, baseline=False,
                 norm_to_file=None, slice_start=None, slice_end=None,
                 energy_min=1e-5, energy_max=1000, energy_step=0.01):
        super().__init__(spectra_file, data_file, repeat, folder)
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
        exp_x_interp = self.exp_x_interp
        exp_y_interp = self.exp_y_interp
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

        # Use lmfit to obtain 'density' to minimize 'y_gap_for_fitting'
        self.fit_result = minimize(y_gap_for_fitting, params_for_fit, method='leastsq',
                                   args=(exp_x_interp, exp_y_interp, self.layer_list,
                                         self.energy_min, self.energy_max, self.energy_step, each_step))
        # Print chi^2
        self.fitted_residual = self.fit_result.__dict__['residual']
        print("Fitting chi^2 : {}".format(sum(self.fitted_residual ** 2)))
        # Print values give best fit
        self.fit_result.__dict__['params'].pretty_print()

        '''Unload fitted params'''
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
        return self.fit_result

    def molar_conc(self):
        self.fitted_simulation = Simulation(energy_min=self.energy_min,
                                            energy_max=self.energy_max,
                                            energy_step=self.energy_step)
        for each_layer in self.layer_list:
            self.fitted_simulation.add_layer(layer=each_layer,
                                             layer_thickness_mm=self.fitted_layer.info[each_layer]['thickness'][
                                                 'value'],
                                             layer_density_gcm3=self.fitted_layer.info[each_layer]['density']['value'])
        for _each_layer in self.layer_list:
            molar_mass_value = self.fitted_simulation.o_reso.stack[_each_layer][_each_layer]['molar_mass']['value']
            molar_mass_units = self.fitted_simulation.o_reso.stack[_each_layer][_each_layer]['molar_mass']['units']
            # Adding molar_mass to fitted_layer info
            self.fitted_layer.info[_each_layer]['molar_mass']['value'] = molar_mass_value
            self.fitted_layer.info[_each_layer]['molar_mass']['units'] = molar_mass_units
            # Adding molar_concentration to fitted_layer info
            molar_conc_value = self.fitted_layer.info[_each_layer]['density']['value'] / molar_mass_value
            molar_conc_units = 'mol/cm3'
            self.fitted_layer.info[_each_layer]['molar_conc']['value'] = molar_conc_value
            self.fitted_layer.info[_each_layer]['molar_conc']['units'] = molar_conc_units
            print('Molar conc. of {} is: {} ({})'.format(_each_layer, molar_conc_value, molar_conc_units))

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
        simu_label = 'ideal'
        exp_label = 'exp'
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
                                                 layer_thickness_mm=self.fitted_layer.info[each_layer]['thickness']['value'],
                                                 layer_density_gcm3=self.fitted_layer.info[each_layer]['density']['value'])
        simu_x, simu_y = self.fitted_simulation.xy_simu(angstrom=False, transmission=False)

        # Get plot labels
        simu_label = 'ideal'
        exp_label = 'exp'
        exp_interp_label = 'exp_interp'
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
