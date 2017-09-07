import glob
import os
import numbers
import re
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.constants import Avogadro
from ImagingReso import _utilities
import os


class Simulation(object):
    pass

class Experiment(object):
    folder = ''
    spectra = ''
    data = ''
    delay_us = float
    source_to_detector_m = float
    repeat = int

    def __init__(self, spectra, data, folder='data', repeat=1, delay_us=0, source_to_detector_m=16.12,
                 angstrom=False, transmission=False):
        if os.path.isdir(folder) is False:
            raise ValueError('Folder specified does not exist')
        # Spectra file
        if os.path.exists(folder+'/'+spectra) is False:
            raise ValueError("Can not find spectra file '{}' in '/{}' folder".format(spectra, folder))
        suffix = spectra[-4:]
        if suffix != '.txt': #or suffix != '.csv':
            raise ValueError("Spectra file must be '.txt' or '.csv'")
        # Data file
        if os.path.exists(folder+'/'+data) is False:
            raise ValueError("Can not find data file '{}' in '/{}' folder".format(data, folder))
        suffix = data[-4:]
        if suffix != '.txt': # or suffix != '.csv':
            raise ValueError("Data file must be '.txt' or '.csv'")
        self.folder = folder
        self.repeat = repeat
        self.source_to_detector_m = source_to_detector_m
        self.delay_us = delay_us
        _file_path = os.path.abspath(os.path.dirname(__file__))
        # Spectra file
        self.path_spectra = os.path.join(_file_path, folder + '/' + spectra)
        # Data file
        self.path_data = os.path.join(_file_path, folder + '/' + data)

        spectra = pd.read_csv(self.path_spectra, sep='\t', header=None)
        self.x_s = np.array(spectra[0])
        self.OB_counts = np.array(spectra[1])
        data = pd.read_csv(self.path_data, sep='\t', header=None)
        if np.array(data[0])[:3] == [1, 2, 3, 4]:
            z = np.array(data[1]) / repeat
        else:
            z = np.array(data[0]) / repeat

        # Output y of experimental data
        if transmission is True:
            self.y = z
        else:
            self.y = 1-z
        # Output x of experimental data
        if angstrom is True:
            self.x = _utilities.s_to_angstroms(self.x_s,
                                               delay_us=delay_us,
                                               source_to_detector_m=source_to_detector_m)
        else:
            self.x = _utilities.s_to_ev(self.x_s,
                                        delay_us=delay_us,
                                        source_to_detector_m=source_to_detector_m)

    # def x_ev(self):
    #     s = self.x_s
    #     delay_us = self.delay_us
    #     source_to_detector_m = self.source_to_detector_m
    #     x_ev = _utilities.s_to_ev(s,delay_us=delay_us,source_to_detector_m=source_to_detector_m)
    #     return x_ev
    #
    # def x_angstrom(self):
    #     s = self.x_s
    #     delay_us = self.delay_us
    #     source_to_detector_m = self.source_to_detector_m
    #     x_angstrom = _utilities.s_to_angstroms(s,delay_us=delay_us,source_to_detector_m=source_to_detector_m)
    #     return x_angstrom





# def get_spectra(_filename, delay_us, source_to_detector_cm, time_lamda_ev_axis='eV'):
#     df_spectra = pd.read_csv(_filename, sep='\t', header=None)
#     time_array = (np.array(df_spectra[0]))
#     # flux_array = (np.array(df_spectra[1]))
#     if time_lamda_ev_axis == 'lamda':
#         lamda_array = _utilities.s_to_angstroms(time_array, delay_us, source_to_detector_cm)
#         return lamda_array
#     if time_lamda_ev_axis == 'eV':
#         ev_array = _utilities.s_to_ev(time_array, delay_us, source_to_detector_cm)
#         return ev_array
#     if time_lamda_ev_axis == 'time':
#         return time_array
#
# def get_spectra2(_filename, delay_us, source_to_detector_cm, time_lamda_ev_axis='eV'):
#     df_spectra = pd.read_csv(_filename, sep='\t', header=None)
#     time_array = (np.array(df_spectra[0]))
#     counts_array = (np.array(df_spectra[1]))
#     if time_lamda_ev_axis == 'lamda':
#         lamda_array = _utilities.s_to_angstroms(time_array, delay_us, source_to_detector_cm)
#         return lamda_array, counts_array
#     if time_lamda_ev_axis == 'eV':
#         ev_array = _utilities.s_to_ev(time_array, delay_us, source_to_detector_cm)
#         return ev_array, counts_array
#     if time_lamda_ev_axis == 'time':
#         return time_array, counts_array
#
#
# def get_spectra_slice(_filename, time_lamda_ev_axis, delay_us, source_to_detector_cm, _slice):
#     df_spectra = pd.read_csv(_filename, sep='\t', header=None)
#     time_array = (np.array(df_spectra[0]))
#     # flux_array = (np.array(df_spectra[1]))
#     if time_lamda_ev_axis == 'lamda':
#         lamda_array = _utilities.s_to_angstroms(time_array, delay_us, source_to_detector_cm)
#         return lamda_array
#     if time_lamda_ev_axis == 'eV':
#         ev_array = _utilities.s_to_ev(time_array, delay_us, source_to_detector_cm)
#         ev_array = ev_array[_slice:]
#         return ev_array
#     if time_lamda_ev_axis == 'lamda':
#         return time_array
