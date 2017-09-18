import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
from scipy.interpolate import interp1d
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation
import numpy as np
from lmfit import minimize


class Calibration(Simulation):
    def __init__(self, spectra, data, layer_1, thickness_1, energy_min, energy_max, energy_step=0.01,
                 repeat=1, density_1=np.NaN, folder='data'):
        super().__init__(layer_1, thickness_1, density_1, energy_min, energy_max, energy_step)
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.energy_step = energy_step
        self.experiment = Experiment(spectra=spectra, data=data, repeat=repeat, folder=folder)
        self.repeat = repeat

    def cost(self, params_exp):
        x_simu, y_simu = self.x, self.y
        source_to_detector_m = params_exp['source_to_detector_m']
        offset_us = params_exp['offset_us']
        x_exp, y_exp = self.experiment.xy_scaled(energy_min=self.energy_min,
                                                 energy_max=self.energy_max,
                                                 energy_step=self.energy_step,
                                                 angstrom=False,
                                                 transmission=False,
                                                 offset_us=offset_us,
                                                 source_to_detector_m=source_to_detector_m)
        # print(x_simu - x_exp)
        # if x_exp != x_simu:
        #     raise ValueError('Energy range and/or energy step entered need to be identical for both simulation and experiment.')

        chi = y_exp - y_simu
        return sum(chi ** 2)

    # def exp_params(self, params_exp):
    #     out = minimize(self.cost(params_exp), params_exp, method='leastsq')
    #
    #     return out
