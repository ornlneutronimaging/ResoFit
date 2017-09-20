import re
import numpy as np
import pandas as pd
import ImagingReso._utilities as reso_utils
from ImagingReso.resonance import Resonance
import os
from lmfit import Parameters
from scipy.interpolate import interp1d
import peakutils as pku
import numbers


class Experiment(object):
    source_to_detector_m = np.NaN
    offset_us = np.NaN
    repeat = np.NaN
    spectra_path = ''
    data_path = ''
    spectra = pd.DataFrame
    data = pd.DataFrame

    def __init__(self, spectra_file, data_file, repeat=1, folder='data'):

        _file_path = os.path.abspath(os.path.dirname(__file__))
        _folder_path = os.path.join(_file_path, folder)

        # Check folder existence
        if os.path.isdir(_folder_path) is False:
            raise ValueError("Folder '{}' specified does not exist".format(folder))
        # Check Spectra file format and existence
        spectra_format = spectra_file[-4:]
        if spectra_format not in ['.txt', '.csv']:
            raise ValueError("Spectra file must be in the format of '.txt' or '.csv'")
        self.spectra_path = os.path.join(_folder_path, spectra_file)
        if os.path.exists(self.spectra_path) is False:
            raise ValueError("Can not find spectra file '{}' in '/{}' folder".format(spectra_file, folder))
        # Data file format and existence
        data_format = data_file[-4:]
        if data_format not in ['.txt', '.csv']:
            raise ValueError("Spectra file must be in the format of '.txt' or '.csv'")
        self.data_path = os.path.join(_folder_path, data_file)
        if os.path.exists(self.data_path) is False:
            raise ValueError("Can not find data file '{}' in '/{}' folder".format(data_file, folder))
        # Force repeat to be an int >=1
        if isinstance(repeat, int) is False:
            raise ValueError("Repeat value must be an integer!")
        if repeat < 1:
            raise ValueError("Repeat value must be an integer >= 1 !")

        # Load spectra file
        self.spectra = pd.read_csv(self.spectra_path, sep='\t', header=None)
        if len(self.spectra.columns) > 2:
            raise ValueError("Spectra file must only contains '2' columns (x & y) !")
        if len(self.spectra.columns) == 1:
            self.spectra = pd.read_csv(self.spectra_path, sep=',', header=None)
        if isinstance(self.spectra[0][0], numbers.Number) is False:
            self.spectra = pd.read_csv(self.data_path, sep='\t', header=None, skiprows=1)
        if len(self.spectra.columns) == 1:
            self.spectra = pd.read_csv(self.spectra_path, sep=',', header=None)
        if isinstance(self.spectra[0][0], numbers.Number) is False:
            raise ValueError("Spectra file contains more than '1' row of descriptions !")

        # Load data file
        self.data = pd.read_csv(self.data_path, sep='\t', header=None)
        if len(self.data.columns) > 2:
            raise ValueError("Data file must only contains '2' columns (x & y) !")
        if len(self.data.columns) == 1:
            self.data = pd.read_csv(self.data_path, sep=',', header=None)
        if isinstance(self.data[0][0], numbers.Number) is False:
            self.data = pd.read_csv(self.data_path, sep='\t', header=None, skiprows=1)
        if len(self.data.columns) == 1:
            self.data = pd.read_csv(self.data_path, sep=',', header=None)
        if isinstance(self.data[0][0], numbers.Number) is False:
            raise ValueError("Data file contains more than '1' row of descriptions !")

    def x_raw(self, angstrom=False, offset_us=0., source_to_detector_m=15.):
        self.offset_us = offset_us
        self.source_to_detector_m = source_to_detector_m
        x_exp_raw = reso_utils.s_to_ev(self.spectra[0],  # x in seconds
                                       offset_us=offset_us,
                                       source_to_detector_m=source_to_detector_m)
        if angstrom is True:
            x_exp_raw = reso_utils.ev_to_angstroms(x_exp_raw)
        return x_exp_raw

    def y_raw(self, transmission=False):
        if len(self.data.columns) == 1:
            y_exp_raw = np.array(self.data[0]) / self.repeat
        else:
            y_exp_raw = np.array(self.data[1]) / self.repeat
        if transmission is False:
            y_exp_raw = 1 - y_exp_raw
        return y_exp_raw

    def xy_scaled(self, energy_min, energy_max, energy_step, angstrom=False, transmission=False,
                  offset_us=0, source_to_detector_m=15.):
        self.offset_us = offset_us
        self.source_to_detector_m = source_to_detector_m
        x_exp_raw = reso_utils.s_to_ev(self.spectra[0],  # x in seconds
                                       offset_us=offset_us,
                                       source_to_detector_m=source_to_detector_m)
        if len(self.data.columns) == 1:
            y_exp_raw = np.array(self.data[0]) / self.repeat
        else:
            y_exp_raw = np.array(self.data[1]) / self.repeat
        if transmission is False:
            y_exp_raw = 1 - y_exp_raw

        nbr_point = (energy_max - energy_min) / energy_step
        x_interp = np.linspace(energy_min, energy_max, nbr_point)
        y_interp_function = interp1d(x=x_exp_raw, y=y_exp_raw, kind='cubic')
        y_interp = y_interp_function(x_interp)
        # baseline = pku.baseline(y_interp)
        # y_interp = y_interp - baseline

        if angstrom is True:
            x_interp = reso_utils.ev_to_angstroms(x_interp)
        return x_interp, y_interp

