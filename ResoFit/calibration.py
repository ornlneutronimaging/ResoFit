import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
from scipy.interpolate import interp1d
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation
import numpy as np
from lmfit import minimize
from ResoFit._gap_functions import y_gap_for_calibration


class Calibration(Simulation):
    def __init__(self, spectra_file, data_file, raw_layer,
                 energy_min=1e-5, energy_max=1000, energy_step=0.01,
                 repeat=1, folder='data', baseline=False):
        """
        Initialization with passed file location and sample info
        :param spectra_file:
        :param data_file:
        :param layer:
        :param thickness_mm:
        :param density_gcm3:
        :param energy_min:
        :param energy_max:
        :param energy_step:
        :param repeat:
        :param folder:
        :param baseline: boolean. True -> to remove baseline/background by detrend
        """
        super().__init__(energy_min=energy_min, energy_max=energy_max, energy_step=energy_step)
        for _each_layer in list(raw_layer.info.keys()):
            self.add_layer(layer=_each_layer,
                           layer_thickness_mm=raw_layer.info[_each_layer]['thickness']['value'],
                           layer_density_gcm3=raw_layer.info[_each_layer]['density']['value'])
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.energy_step = energy_step
        self.experiment = Experiment(spectra_file=spectra_file, data_file=data_file, repeat=repeat, folder=folder)
        self.repeat = repeat
        self.data_file = data_file
        self.spectra_file = spectra_file
        self.init_source_to_detector_m = None
        self.init_offset_us = None
        self.calibrated_offset_us = None
        self.calibrated_source_to_detector_m = None
        self.calibrate_result = None
        self.exp_x_raw_calibrated = None
        self.exp_y_raw_calibrated = None
        self.exp_x_interp_calibrated = None
        self.exp_y_interp_calibrated = None
        self.baseline = baseline
        self.calibrated_residual = None
        self.raw_layer = raw_layer

    def norm_to(self, file):
        if file is not None:
            self.experiment.norm_to(file=file)

    def slice(self, slice_start=None, slice_end=None):
        self.experiment.slice(slice_start=slice_start, slice_end=slice_end)

    def calibrate(self, source_to_detector_m, offset_us, vary='all', each_step=False):
        """
        calibrate the instrumental parameters: source-to-detector-distance & detector delay
        :param each_step: boolean. True -> show values and chi^2 of each step
        :param source_to_detector_m: estimated distance in m
        :param offset_us: estimated time offset in us
        :param vary: vary one of or both of 'source_to_detector' and 'offset' to calibrate (default: 'all')

        :return: lmfit MinimizerResult
        """
        self.init_source_to_detector_m = source_to_detector_m
        self.init_offset_us = offset_us
        if vary not in ['source_to_detector', 'offset', 'all', 'none']:
            raise ValueError("'vary=' can only be one of ['source_to_detector', 'offset', 'all' 'none']")
        simu_x = self.simu_x
        simu_y = self.simu_y

        source_to_detector_vary_tag = True
        offset_vary_tag = True
        if vary == 'source_to_detector':
            offset_vary_tag = False
        if vary == 'offset':
            source_to_detector_vary_tag = False
        if vary == 'none':
            source_to_detector_vary_tag = False
            offset_vary_tag = False
        params_calibrate = Parameters()
        params_calibrate.add('source_to_detector_m', value=source_to_detector_m, vary=source_to_detector_vary_tag)
        params_calibrate.add('offset_us', value=offset_us, vary=offset_vary_tag)
        # Print before
        print("Params before calibration:")
        params_calibrate.pretty_print()
        # Use lmfit to obtain 'source_to_detector_m' & 'offset_us' to minimize 'y_gap_for_calibration'
        self.calibrate_result = minimize(y_gap_for_calibration, params_calibrate, method='leastsq',
                                         args=(simu_x, simu_y,
                                               self.energy_min, self.energy_max, self.energy_step,
                                               self.experiment, self.baseline, each_step))
        # Print after
        print("Params after calibration:")
        self.calibrate_result.__dict__['params'].pretty_print()
        # Print chi^2
        self.calibrated_residual = self.calibrate_result.__dict__['residual']
        print("Calibration chi^2 : {}\n".format(sum(self.calibrated_residual ** 2)))
        self.calibrated_offset_us = self.calibrate_result.__dict__['params'].valuesdict()['offset_us']
        self.calibrated_source_to_detector_m = \
            self.calibrate_result.__dict__['params'].valuesdict()['source_to_detector_m']

        # Save the calibrated experimental x & y in Calibration class
        self.exp_x_raw_calibrated = self.experiment.x_raw(angstrom=False,
                                                          offset_us=self.calibrated_offset_us,
                                                          source_to_detector_m=self.calibrated_source_to_detector_m)
        self.exp_y_raw_calibrated = self.experiment.y_raw(transmission=False, baseline=self.baseline)

        self.exp_x_interp_calibrated, self.exp_y_interp_calibrated = self.experiment.xy_scaled(
            energy_min=self.energy_min,
            energy_max=self.energy_max,
            energy_step=self.energy_step,
            offset_us=self.calibrated_offset_us,
            source_to_detector_m=self.calibrated_source_to_detector_m,
            baseline=self.baseline)

        return self.calibrate_result

    def plot(self, interp=False, before=False):
        """
        Plot the raw experimental data and theoretical resonance signal after calibration
        :param before: boolean. Plot the data before calibration applied.
        :param interp: boolean. True -> display interpolated exp data
                                False -> display raw exp data
        :return: plot of raw experimental data and theoretical resonance signal after calibration
        """
        simu_label = 'Ideal'
        exp_label = 'Exp'
        exp_before_label = 'Exp_before_calibration'
        exp_interp_label = 'Exp_interp'
        for each_layer in self.layer_list:
            simu_label = simu_label + '_' + each_layer
            exp_label = exp_label + '_' + each_layer
            exp_interp_label = exp_interp_label + '_' + each_layer
            exp_before_label = exp_before_label + '_' + each_layer
        plt.plot(self.simu_x, self.simu_y, 'b-', label=simu_label, markersize=1)
        if interp is False:
            plt.plot(self.exp_x_raw_calibrated, self.exp_y_raw_calibrated, 'ro', label=exp_label, markersize=1)
        else:
            plt.plot(self.exp_x_interp_calibrated, self.exp_y_interp_calibrated, 'r-.', label=exp_interp_label, markersize=1)
        if before is True:
            plt.plot(self.experiment.x_raw(offset_us=self.init_offset_us,
                                           source_to_detector_m=self.init_source_to_detector_m),
                     self.experiment.y_raw(baseline=self.baseline),
                     'ko', label=exp_before_label, markersize=1, alpha=0.6)

        plt.title('Calibration result')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Attenuation')
        plt.ylim(ymax=1.01)
        plt.xlim(0, self.energy_max)
        plt.legend(loc='best')
        plt.show()
