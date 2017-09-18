import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
from scipy.interpolate import interp1d
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation
import numpy as np
from lmfit import minimize
from ResoFit._utilities import peak_y_gap


class Calibration(Simulation):
    def __init__(self, spectra, data, layer_1, thickness_1, density_1=np.NaN,
                 energy_min=1e-5, energy_max=1000, energy_step=0.01,
                 repeat=1, folder='data'):
        super().__init__(layer_1, thickness_1, density_1, energy_min, energy_max, energy_step)
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.energy_step = energy_step
        self.experiment = Experiment(spectra=spectra, data=data, repeat=repeat, folder=folder)
        self.repeat = repeat
        self.data = data
        self.spectra = spectra

    def get_exp_params(self, params_init):
        simu_x = self.x
        simu_y = self.y,
        energy_min = self.energy_min
        energy_max = self.energy_max
        energy_step = self.energy_step
        data = self.data
        spectra = self.spectra
        repeat = self.repeat
        out = minimize(peak_y_gap(params_init, simu_x, simu_y, energy_min, energy_max, energy_step, data, spectra, repeat),
                       params_init, method='leastsq',
                       args=(simu_x, simu_y, energy_min, energy_max, energy_step, data, spectra, repeat))
        print(out)
