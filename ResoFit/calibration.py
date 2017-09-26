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
    def __init__(self, spectra_file, data_file, layer_1, thickness_1, density_1=np.NaN,
                 energy_min=1e-5, energy_max=1000, energy_step=0.01,
                 repeat=1, folder='data'):
        """
        Initialization with passed file location and sample info
        :param spectra_file:
        :param data_file:
        :param layer_1:
        :param thickness_1:
        :param density_1:
        :param energy_min:
        :param energy_max:
        :param energy_step:
        :param repeat:
        :param folder:
        """
        super().__init__(energy_min, energy_max, energy_step)
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.energy_step = energy_step
        self.add_layer(layer=layer_1, layer_thickness=thickness_1, layer_density=density_1)
        self.experiment = Experiment(spectra_file=spectra_file, data_file=data_file, repeat=repeat, folder=folder)
        self.repeat = repeat
        self.data_file = data_file
        self.spectra_file = spectra_file
        self.calibrated_offset_us = None
        self.calibrated_source_to_detector_m = None
        self.calibrate_result = None
        self.exp_x_raw_calibrated = None
        self.exp_y_raw_calibrated = None
        self.exp_x_interp_calibrated = None
        self.exp_y_interp_calibrated = None
        self.layer_1 = layer_1

    def calibrate(self, source_to_detector_m, offset_us, vary='all'):
        """
        calibrate the instrumental parameters: source-to-detector-distance & detector delay
        :param source_to_detector_m: estimated distance in m
        :param offset_us: estimated time offset in us
        :param vary: vary one of or both of 'source_to_detector' and 'offset' to calibrate (default: 'all')
        :return: lmfit MinimizerResult
        """
        if vary not in ['source_to_detector', 'offset', 'all']:
            raise ValueError("'vary=' can only be one of ['source_to_detector', 'offset', 'all']")
        simu_x = self.simu_x
        simu_y = self.simu_y

        source_to_detector_vary_tag = True
        offset_vary_tag = True
        if vary == 'source_to_detector':
            offset_vary_tag = False
        if vary == 'offset':
            source_to_detector_vary_tag = False
        params_calibrate = Parameters()
        params_calibrate.add('source_to_detector_m', value=source_to_detector_m, vary=source_to_detector_vary_tag)
        params_calibrate.add('offset_us', value=offset_us, vary=offset_vary_tag)
        # Use lmfit to obtain 'source_to_detector_m' & 'offset_us' to minimize 'y_gap_for_calibration'
        self.calibrate_result = minimize(y_gap_for_calibration, params_calibrate, method='leastsq',
                                         args=(simu_x, simu_y,
                                               self.energy_min, self.energy_max, self.energy_step,
                                               self.data_file, self.spectra_file, self.repeat))
        self.calibrated_offset_us = self.calibrate_result.__dict__['params'].valuesdict()['offset_us']
        self.calibrated_source_to_detector_m = \
            self.calibrate_result.__dict__['params'].valuesdict()['source_to_detector_m']
        # Print values give best fit
        self.calibrate_result.__dict__['params'].pretty_print()
        # Save the calibrated experimental x & y in Calibration class
        self.exp_x_raw_calibrated = self.experiment.x_raw(angstrom=False,
                                                          offset_us=self.calibrated_offset_us,
                                                          source_to_detector_m=self.calibrated_source_to_detector_m)
        self.exp_y_raw_calibrated = self.experiment.y_raw(transmission=False)

        self.exp_x_interp_calibrated, self.exp_y_interp_calibrated = self.experiment.xy_scaled(
            energy_min=self.energy_min,
            energy_max=self.energy_max,
            energy_step=self.energy_step,
            offset_us=self.calibrated_offset_us,
            source_to_detector_m=self.calibrated_source_to_detector_m)

        return self.calibrate_result

    def plot_before(self):
        """
        Plot the raw experimental data and theoretical resonance signal before calibration
        :return:
        """
        plt.plot(self.simu_x, self.simu_y,
                 'b.', label=self.layer_1 + '_ideal', markersize=1)

        plt.plot(self.experiment.x_raw(), self.experiment.y_raw(),
                 'r.', label=self.layer_1 + '_exp', markersize=1)

        plt.title('Before Calibration')
        plt.ylim(-0.01, 1.01)
        plt.xlim(0, self.energy_max)
        plt.legend(loc='best')
        plt.show()

    def plot_after(self):
        """
        Plot the raw experimental data and theoretical resonance signal after calibration
        :return:
        """
        plt.plot(self.simu_x, self.simu_y,
                 'b.', label=self.layer_1 + '_ideal', markersize=1)

        plt.plot(self.exp_x_raw_calibrated, self.exp_y_raw_calibrated,
                 'r.', label=self.layer_1 + '_exp', markersize=1)

        plt.title('After Calibration')
        plt.ylim(-0.01, 1.01)
        plt.xlim(0, self.energy_max)
        plt.legend(loc='best')
        plt.show()

    def plot_after_interp(self):
        """
        Plot the interpolated experimental data and theoretical resonance signal after calibration
        :return:
        """
        plt.plot(self.simu_x, self.simu_y,
                 'b.', label=self.layer_1 + '_ideal', markersize=1)

        plt.plot(self.exp_x_interp_calibrated, self.exp_y_interp_calibrated,
                 'r.', label=self.layer_1 + '_exp', markersize=1)

        plt.title('After Calibration')
        plt.ylim(-0.01, 1.01)
        plt.xlim(0, self.energy_max)
        plt.legend(loc='best')
        plt.show()
