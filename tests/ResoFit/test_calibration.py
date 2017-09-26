import unittest
from ResoFit.calibration import Calibration
import numpy as np
import os


class TestCalibration(unittest.TestCase):
    def setUp(self):
        _file_path = os.path.dirname(__file__)
        self.folder = os.path.abspath(os.path.join(_file_path, '../../ResoFit/data/_mock_data_for_test'))
        self.data_file = os.path.join(self.folder, '_data_calibration_fitting_test.txt')
        self.spectra_file = os.path.join(self.folder, '_spectra_calibration_fitting_test.txt')

        energy_min = 7
        energy_max = 150
        energy_step = 0.01

        thickness = 0.075
        layer = 'Gd'
        self.calibration = Calibration(data_file=self.data_file,
                                       spectra_file=self.spectra_file,
                                       layer_1=layer,
                                       thickness_1=thickness,
                                       density_1=np.NaN,
                                       energy_min=energy_min,
                                       energy_max=energy_max,
                                       energy_step=energy_step,
                                       folder=self.folder)

    source_to_detector_m = 16.
    offset_us = 0

    def test_vary_tag(self):
        """assert vary to be one of ['source_to_detector', 'offset', 'all']"""
        data_file = self.data_file
        spectra_file = self.spectra_file
        self.assertRaises(ValueError, self.calibration.calibrate,
                          source_to_detector_m=self.source_to_detector_m,
                          offset_us=self.offset_us,
                          vary='tag_not_exist')
