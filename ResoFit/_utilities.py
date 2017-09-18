import re
import numpy as np
import pandas as pd
from ImagingReso import _utilities
from ImagingReso.resonance import Resonance
import os
from lmfit import Parameters
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation


def y_gap_for_calibration(params, simu_x, simu_y, energy_min, energy_max, energy_step, data_file, spectra_file,
                          repeat=1):
    # Unpack Parameters:
    parvals = params.valuesdict()
    source_to_detector_m = parvals['source_to_detector_m']
    offset_us = parvals['offset_us']
    experiment = Experiment(data_file=data_file, spectra_file=spectra_file, repeat=repeat)
    exp_x, exp_y = experiment.xy_scaled(energy_min=energy_min,
                                        energy_max=energy_max,
                                        energy_step=energy_step,
                                        angstrom=False,
                                        transmission=False,
                                        offset_us=offset_us,
                                        source_to_detector_m=source_to_detector_m)
    # if sum((simu_x - exp_x) ** 2) >= 0.001:
    #     raise ValueError("The experiment x-axis is not identical to simulation x-axis!")
    gap = (exp_y - simu_y) ** 2
    return gap


def y_gap_for_fitting(params, exp_x, exp_y, layer, energy_min, energy_max, energy_step):
    parvals = params.valuesdict()
    density = parvals['density']
    thickness = parvals['thickness']
    simulation = Simulation(layer_1=layer,
                            density_1=density,
                            thickness_1=thickness,
                            energy_min=energy_min,
                            energy_max=energy_max,
                            energy_step=energy_step)
    pass
