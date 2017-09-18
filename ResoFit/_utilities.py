import re
import numpy as np
import pandas as pd
from ImagingReso import _utilities
from ImagingReso.resonance import Resonance
import os
from lmfit import Parameters
from ResoFit.experiment import Experiment


def cost(simu_x, simu_y, params):
    # source_to_detector_m = params['source_to_detector_m']
    # offset_us = params['offset_us']
    # experiment = Experiment(data='all_thin.txt', spectra='Image002_Spectra.txt', repeat=5, offset_us=offset_us, source_to_detector_m=source_to_detector_m)
    # exp_x = experiment.x
    # baseline = pku.baseline(experiment.y)
    # exp_y = experiment.y - baseline
    # exp_y_function = interp1d(x=exp_x, y=exp_y, kind='cubic')
    # exp_y_interp = exp_y_function(simu_x)
    #
    # simulation = Simulation(layer_1=_layer_1,
    #                         thickness_1=_thickness_1,
    #                         density_1=np.NaN,
    #                         _energy_min=_energy_min,
    #                         _energy_max=_energy_max,
    #                         _energy_step=_energy_step)
    # simu_x = simulation.x()
    # simu_y = simulation.y()
    #
    #
    # chi = exp_y_interp - simu_y
    # return chi**2
    pass


# def peak_x_gap(params, ideal_x_index, y_data_array):
#     # Unpack Parameters:
#     parvals = params.valuesdict()
#     source_to_detector_cm = parvals['source_to_detector_cm']
#     delay_us = parvals['delay_us']
#     # Model:
#     spectra_path = 'data/spectra.txt'
#     range_min = 500
#     range_max = 2000
#     x_data_array = _functions.get_spectra_range(spectra_path, delay_us,
#                                                 source_to_detector_cm, range_min, range_max)
#     exp_y_index = pku.indexes(y_data_array, thres=0.12/max(y_data_array), min_dist=7)
#     exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
#     # gap = (exp_x_index[0] - ideal_x_index) ** 2
#     # print(exp_x_index)
#     # print(ideal_x_index)
#     gap = (exp_x_index - ideal_x_index) ** 2
#     return gap


def peak_y_gap(params, simu_x, simu_y, energy_min, energy_max, energy_step, data_file, spectra_file, repeat=1):
    # Unpack Parameters:
    parvals = params.valuesdict()
    source_to_detector_m = parvals['source_to_detector_m']
    offset_us = parvals['offset_us']
    experiment = Experiment(data=data_file, spectra=spectra_file, repeat=repeat)
    exp_x, exp_y = experiment.xy_scaled(energy_min=energy_min,
                                        energy_max=energy_max,
                                        energy_step=energy_step,
                                        angstrom=False,
                                        transmission=False,
                                        offset_us=offset_us,
                                        source_to_detector_m=source_to_detector_m)
    print(simu_x is exp_x)
    gap = (exp_y - simu_y)**2
    return gap

    # Model:


    # def peak_x_gap_scipy(delay_us, ideal_x_index, y_data_array):
    #     # Unpack Parameters:
    #     source_to_detector_cm = 1610.9  # cm
    #     # Model:
    #     time_lamda_ev_axis = 'eV'
    #     spectra_path = 'data/spectra.txt'
    #     _slice = 220
    #     x_data_array = _functions.get_spectra_slice(spectra_path, time_lamda_ev_axis, delay_us,
    #                                                 source_to_detector_cm, _slice)
    #     exp_y_index = pku.indexes(y_data_array, thres=0.6, min_dist=50)
    #     exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
    #     gap = (exp_x_index[0] - ideal_x_index) ** 2
    #     return gap


    # def peak_gap(params, ideal_x_index):
    #     # Unpack Parameters:
    #     parvals = params.valuesdict()
    #     source_to_detector_cm = parvals['source_to_detector_cm']
    #     time_us = parvals['time_us']
    #     # Model:
    #     energy_miliev = 81.787 / (0.3956 * time_us / source_to_detector_cm) ** 2
    #     energy_ev = energy_miliev / 1000
    #     return (energy_ev - ideal_x_index) ** 2
