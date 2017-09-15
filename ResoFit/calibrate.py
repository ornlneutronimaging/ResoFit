import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
from scipy.interpolate import interp1d
from ResoFit.experiment import Experiment


class Calibrate(Experiment):

    def offset_us(self, energy_min=1e-5, energy_max=1000, energy_step=0.01):

        self.params_exp = Parameters()
        self.params_exp.add('source_to_detector_m', value=self.source_to_detector_m)
        self.params_exp.add('offset_us', value=self.offset_us)

        _source_to_detector_m = self.params_exp['source_to_detector_m']
        _offset_us = self.params_exp['offset_us']
        experiment = Experiment(data='all_thin.txt', spectra='Image002_Spectra.txt', repeat=5,
                                source_to_detector_m=_source_to_detector_m, offset_us=_offset_us)
        exp_x, exp_y = experiment.xy_scaled(energy_min, energy_max, energy_step)
        baseline = pku.baseline(exp_y)
        exp_y = exp_y - baseline
        chi = exp_y_interp - simu_
        return chi ** 2
