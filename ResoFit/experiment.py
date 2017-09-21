import numpy as np
import pandas as pd
import ImagingReso._utilities as reso_utils
import os
from scipy.interpolate import interp1d
import peakutils as pku
from ResoFit._utilities import load_txt_csv


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
        # Check the folder existence
        if os.path.isdir(_folder_path) is False:
            raise ValueError("Folder '{}' specified does not exist".format(folder))
        # Spectra file format and existence
        spectra_format = spectra_file[-4:]
        if spectra_format not in ['.txt', '.csv']:
            raise ValueError("Spectra file must be in the format of '.txt' or '.csv'")
        self.spectra_path = os.path.join(_folder_path, spectra_file)
        if os.path.exists(self.spectra_path) is False:
            raise ValueError("Can not find spectra file '{}' in '{}' folder".format(spectra_file, _folder_path))
        # Data file format and existence
        date_format = data_file[-4:]
        if date_format not in ['.txt', '.csv']:
            raise ValueError("Spectra file must be in the format of '.txt' or '.csv'")
        self.data_path = os.path.join(_folder_path, data_file)
        if os.path.exists(self.data_path) is False:
            raise ValueError("Can not find data file '{}' in '{}' folder".format(data_file, _folder_path))
        # Force repeat be an int >=1
        if isinstance(repeat, int) is False:
            raise ValueError("Repeat value must be an integer!")
        if repeat < 1:
            raise ValueError("Repeat value must be an integer >= 1 !")

        # Load data
        self.repeat = repeat
        self.spectra = load_txt_csv(self.spectra_path)
        self.data = load_txt_csv(self.data_path)

        if self.spectra[0][0] is str:
            if self.spectra[0][0].islower() or self.spectra[0][0].isupper() is True:
                raise ValueError("Remove the axis descriptions in '{}' before loading ".format(spectra_file))
            else:
                raise ValueError("The file '{}' columns must be separated with '\t' or ',' ".format(spectra_file))

        if self.data[0][0] is str:
            if self.data[0][0].islower() or self.data[0][0].isupper() is True:
                raise ValueError("Remove the axis descriptions in '{}' before loading ".format(data_file))
            else:
                raise ValueError("The file '{}' columns must be separated with '\t' or ',' ".format(data_file))

    def x_raw(self, angstrom=False, offset_us=0., source_to_detector_m=15.12):
        self.offset_us = offset_us
        self.source_to_detector_m = source_to_detector_m
        x_exp_raw = reso_utils.s_to_ev(self.spectra[0],  # x in seconds
                                       offset_us=offset_us,
                                       source_to_detector_m=source_to_detector_m)
        if angstrom is True:
            x_exp_raw = reso_utils.ev_to_angstroms(x_exp_raw)
        return x_exp_raw

    def y_raw(self, transmission=False):
        if np.array(self.data[0])[:3] == [1, 2, 3, 4]:
            y_exp_raw = np.array(self.data[1]) / self.repeat
        else:
            y_exp_raw = np.array(self.data[0]) / self.repeat
        if transmission is False:
            y_exp_raw = 1 - y_exp_raw
        return y_exp_raw

    def xy_scaled(self, energy_min, energy_max, energy_step, angstrom=False, transmission=False,
                  offset_us=0, source_to_detector_m=15.12):
        self.offset_us = offset_us
        self.source_to_detector_m = source_to_detector_m
        x_exp_raw = reso_utils.s_to_ev(self.spectra[0],  # x in seconds
                                       offset_us=offset_us,
                                       source_to_detector_m=source_to_detector_m)
        if np.array(self.data[0])[:3] == [1, 2, 3, 4]:
            y_exp_raw = np.array(self.data[1]) / self.repeat
        else:
            y_exp_raw = np.array(self.data[0]) / self.repeat
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
