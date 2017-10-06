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


class FitResonance(Experiment):
    fit_result = None
    fitted_density_gcm3 = None
    fitted_thickness_mm = None
    fitted_residual = None
    fitted_gap = None
    fitted_fjac = None

    def __init__(self, spectra_file, data_file,
                 calibrated_offset_us, calibrated_source_to_detector_m,
                 layer, layer_thickness_mm=0.2, layer_density_gcm3=np.NaN,
                 folder='data', repeat=1, baseline=False,
                 norm_to_file=None, slice_start=None, slice_end=None,
                 energy_min=1e-5, energy_max=1000, energy_step=0.01):
        super().__init__(spectra_file, data_file, repeat, folder)
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.energy_step = energy_step
        self.calibrated_offset_us = calibrated_offset_us
        self.calibrated_source_to_detector_m = calibrated_source_to_detector_m
        self.layer = layer
        self.layer_thickness_mm = layer_thickness_mm
        self.layer_density_gcm3 = layer_density_gcm3
        self.slice(slice_start=slice_start, slice_end=slice_end)
        self.baseline = baseline
        if norm_to_file is not None:
            self.norm_to(norm_to_file)
        # self.add_layer(layer=layer_list[0], layer_thickness=thickness_list[0], density)
        self.exp_x_interp, self.exp_y_interp = self.xy_scaled(energy_min=self.energy_min,
                                                              energy_max=self.energy_max,
                                                              energy_step=self.energy_step,
                                                              angstrom=False, transmission=False,
                                                              offset_us=self.calibrated_offset_us,
                                                              source_to_detector_m=self.calibrated_source_to_detector_m,
                                                              baseline=self.baseline)

    def fit(self, thickness_mm, density_gcm3, vary='density', each_step=False):
        if vary not in ['density', 'thickness', 'all', 'none']:
            raise ValueError("'vary=' can only be one of ['density', 'thickness', 'all', 'none']")
        exp_x_interp = self.exp_x_interp
        exp_y_interp = self.exp_y_interp
        if density_gcm3 is np.NaN:
            density_gcm3 = pt.elements.isotope(self.layer).density
        thickness_vary_tag = False
        density_vary_tag = True
        if vary == 'thickness':
            thickness_vary_tag = True
            density_vary_tag = False
        if vary == 'all':
            thickness_vary_tag = True
        if vary == 'none':
            density_vary_tag = False
        params_fit = Parameters()
        params_fit.add('thickness_mm', value=thickness_mm, vary=thickness_vary_tag, min=0)
        params_fit.add('density_gcm3', value=density_gcm3, vary=density_vary_tag, min=0)

        # Use lmfit to obtain 'density' to minimize 'y_gap_for_fitting'
        self.fit_result = minimize(y_gap_for_fitting, params_fit, method='leastsq',
                                   args=(exp_x_interp, exp_y_interp, self.layer,
                                         self.energy_min, self.energy_max, self.energy_step, each_step))
        # Print chi^2
        self.fitted_residual = self.fit_result.__dict__['residual']
        print("Fitting chi^2 : {}".format(sum(self.fitted_residual**2)))
        # Print values give best fit
        self.fit_result.__dict__['params'].pretty_print()
        # Save the fitted 'density' or 'thickness' in FitResonance class
        self.fitted_density_gcm3 = self.fit_result.__dict__['params'].valuesdict()['density_gcm3']
        self.fitted_thickness_mm = self.fit_result.__dict__['params'].valuesdict()['thickness_mm']
        # self.fitted_fjac = self.fit_result.__dict__['fjac']
        # print(self.fit_result.__dict__['fjac'][0])

        return self.fit_result

    def molar_conc(self, element):
        layer = self.layer
        # Check if element exist
        _formula = re.findall(r'([A-Z][a-z]*)(\d*)', layer)
        _elements = []
        for _element in _formula:
            _single_element = list(_element)[0]
            _elements.append(_single_element)
        if element not in _elements:
            raise ValueError('Element {} specified does not exist in {} layer.'.format(element, layer))
        # convert fitted elemental density to molar concentration
        simulation = Simulation(energy_min=self.energy_min,
                                energy_max=self.energy_max,
                                energy_step=self.energy_step)
        simulation.add_layer(layer=self.layer, layer_thickness_mm=self.layer_thickness_mm, layer_density_gcm3=self.layer_density_gcm3)
        molar_mass = simulation.o_reso.stack[layer][element]['molar_mass']['value']
        molar_conc = self.fitted_density_gcm3 / molar_mass
        print('Molar conc. of element {} in layer {} is: {} (mol/cm3)'.format(element, layer, molar_conc))
        return molar_conc

    def plot_before(self):
        simulation = Simulation(energy_min=self.energy_min,
                                energy_max=self.energy_max,
                                energy_step=self.energy_step)
        simulation.add_layer(layer=self.layer, layer_thickness_mm=self.layer_thickness_mm, layer_density_gcm3=self.layer_density_gcm3)
        simu_x, simu_y = simulation.xy_simu(angstrom=False, transmission=False)
        plt.plot(simu_x, simu_y,
                 'b-', label=self.layer + '_simu', markersize=1)

        plt.plot(self.x_raw(angstrom=False, offset_us=self.calibrated_offset_us,
                            source_to_detector_m=self.source_to_detector_m),
                 self.y_raw(transmission=False, baseline=self.baseline),
                 'ro', label=self.layer + '_exp', markersize=1)

        plt.title('Before fitting')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Attenuation')
        plt.ylim(ymax=1.01)
        plt.xlim(0, self.energy_max)
        plt.legend(loc='best')
        plt.show()

    def plot_after(self, error=True):
        simulation = Simulation(energy_min=self.energy_min,
                                energy_max=self.energy_max,
                                energy_step=self.energy_step)
        simulation.add_layer(layer=self.layer, layer_thickness_mm=self.fitted_thickness_mm, layer_density_gcm3=self.fitted_density_gcm3)
        simu_x, simu_y = simulation.xy_simu(angstrom=False, transmission=False)
        plt.plot(simu_x, simu_y,
                 'b-', label=self.layer + '_simu', markersize=1)

        plt.plot(self.x_raw(angstrom=False, offset_us=self.calibrated_offset_us,
                            source_to_detector_m=self.source_to_detector_m),
                 self.y_raw(transmission=False, baseline=self.baseline),
                 'ro', label=self.layer + '_exp', markersize=1)

        if error is True:
            # Plot fitting differences
            plt.plot(simu_x, self.fitted_residual-0.2, 'g-', label='Diff.', alpha=0.8)

        plt.title('Best fit')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Attenuation')
        plt.ylim(ymax=1.01)
        plt.xlim(0, self.energy_max)
        plt.legend(loc='best')
        plt.show()
