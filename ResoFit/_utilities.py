import re
import numpy as np
import pandas as pd
from ImagingReso import _utilities
from ImagingReso.resonance import Resonance
import os
from lmfit import Parameters
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation




def cost(simu_x, simu_y, params):
    source_to_detector_m = params['source_to_detector_m']
    offset_us = params['offset_us']
    experiment = Experiment(data='all_thin.txt', spectra='Image002_Spectra.txt', repeat=5, offset_us=offset_us, source_to_detector_m=source_to_detector_m)
    exp_x = experiment.x
    baseline = pku.baseline(experiment.y)
    exp_y = experiment.y - baseline
    exp_y_function = interp1d(x=exp_x, y=exp_y, kind='cubic')
    exp_y_interp = exp_y_function(simu_x)

    simulation = Simulation(layer_1=_layer_1,
                            thickness_1=_thickness_1,
                            density_1=np.NaN,
                            _energy_min=_energy_min,
                            _energy_max=_energy_max,
                            _energy_step=_energy_step)
    simu_x = simulation.x()
    simu_y = simulation.y()


    chi = exp_y_interp - simu_y
    return chi**2
