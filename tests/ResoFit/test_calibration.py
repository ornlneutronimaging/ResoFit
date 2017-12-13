import unittest
from ResoFit.calibration import Calibration
import numpy as np
import os
from ResoFit._utilities import Layer


class TestCalibration(unittest.TestCase):
    database = '_data_for_unittest'

    def setUp(self):
        _file_path = os.path.dirname(__file__)
        self.folder = os.path.abspath(os.path.join(_file_path, '../../ResoFit/data/_mock_data_for_test'))
        self.data_file = os.path.join(self.folder, '_data_calibration_fitting_test.txt')
        self.spectra_file = os.path.join(self.folder, '_spectra_calibration_fitting_test.txt')

        energy_min = 7
        energy_max = 150
        energy_step = 0.01

        layer_1 = 'U'
        thickness_1 = 0.075
        density_1 = None
        layer = Layer()
        layer.add_layer(layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)

        self.calibration = Calibration(data_file=self.data_file,
                                       spectra_file=self.spectra_file,
                                       raw_layer=layer,
                                       energy_min=energy_min,
                                       energy_max=energy_max,
                                       energy_step=energy_step,
                                       folder=self.folder,
                                       database=self.database)

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
