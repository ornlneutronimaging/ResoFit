import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
from scipy.interpolate import interp1d
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation
import numpy as np


class Calibration(Experiment):
    def __init__(self, spectra, data, layer_1, thickness_1, energy_min, energy_max, energy_step, repeat,
                 density_1=np.NaN):
        super().__init__(spectra, data, repeat)
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.energy_step = energy_step
        simulation = Simulation(layer_1=layer_1, thickness_1=thickness_1, density_1=density_1,
                                energy_min=energy_min, energy_max=energy_max, energy_step=energy_step)
        self.x_simu, self.y_simu = simulation.xy_simu()
        self.repeat = repeat

    def cost(self, params_exp):
        source_to_detector_m = params_exp['source_to_detector_m']
        offset_us = params_exp['offset_us']
        x_exp, y_exp = self.xy_scaled(energy_min=self.energy_min,
                                      energy_max=self.energy_max,
                                      energy_step=self.energy_step,
                                      angstrom=False,
                                      transmission=False,
                                      offset_us=offset_us,
                                      source_to_detector_m=source_to_detector_m)
        chi = y_exp - self.y_simu
        return sum(chi ** 2)
