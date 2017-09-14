import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
from scipy.interpolate import interp1d
from ResoFit.experiment import Experiment


class Calibrate(Experiment):

    def offset_us(self):

        self.params_exp = Parameters()
        self.params_exp.add('source_to_detector_m', value=self.source_to_detector_m)
        self.params_exp.add('offset_us', value=self.offset_us)

        _source_to_detector_m = self.params_exp['source_to_detector_m']
        _offset_us = self.params_exp['offset_us']

        experiment = Experiment(data='all_thin.txt', spectra='Image002_Spectra.txt', repeat=5,
                                source_to_detector_m=_source_to_detector_m, offset_us=_offset_us)
        exp_x = experiment.x
        baseline = pku.baseline(experiment.y)
        exp_y = experiment.y - baseline
        exp_y_function = interp1d(x=exp_x, y=exp_y, kind='cubic')
        exp_y_interp = exp_y_function(simu_x)
        chi = exp_y_interp - simu_y
        return chi ** 2

