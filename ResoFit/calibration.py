import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
from scipy.interpolate import interp1d
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation
import numpy as np
from lmfit import minimize
from ResoFit._utilities import y_gap_for_calibration


class Calibration(Simulation):
    def __init__(self, spectra_file, data_file, layer_1, thickness_1, density_1=np.NaN,
                 energy_min=1e-5, energy_max=1000, energy_step=0.01,
                 repeat=1, folder='data'):
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
        self.layer_1 = layer_1

    def calibrate(self, params_calibrate):
        simu_x = self.simu_x
        simu_y = self.simu_y
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

        return self.calibrate_result

    def plot_before(self):
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
        plt.plot(self.simu_x, self.simu_y,
                 'b.', label=self.layer_1 + '_ideal', markersize=1)

        plt.plot(self.exp_x_raw_calibrated, self.exp_y_raw_calibrated,
                 'r.', label=self.layer_1 + '_exp', markersize=1)

        plt.title('After Calibration')
        plt.ylim(-0.01, 1.01)
        plt.xlim(0, self.energy_max)
        plt.legend(loc='best')
        plt.show()

