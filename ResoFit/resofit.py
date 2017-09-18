import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
from scipy.interpolate import interp1d
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation
import numpy as np
from lmfit import minimize
from ResoFit._utilities import y_gap_for_calibration


class FitResonance(Simulation):
    def __init__(self, layer_1, thickness_1):
        super().__init__(layer_1, thickness_1)

    def fit(self):
        pass

    def plot_before(self):
        pass

    def plot_after(self):
        pass
