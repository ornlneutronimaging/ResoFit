import unittest
import os
import numpy as np
from ResoFit.experiment import Experiment


class TestExperiment(unittest.TestCase):
    def setUp(self):
        _file_path = os.path.dirname(__file__)
        self.folder = os.path.abspath(os.path.join(_file_path, '../../ResoFit/data/_mock_data_for_test'))
        self.data_file = os.path.join(self.folder, '_data_unit_test.txt')
        self.spectra_file = os.path.join(self.folder, '_spectra_unit_test.txt')

    energy_min = 7
    energy_max = 8
    energy_step = 0.01

    def test_folder(self):
        """assert given folder existence"""
        folder = 'folder_not_exist'
        data_file = self.data_file
        spectra_file = self.spectra_file
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_spectra_file_format(self):
        """assert given spectra_file file existence and format"""
        folder = self.folder
        data_file = self.data_file
        spectra_file = self.spectra_file + '.pdf'
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_data_file_format(self):
        """assert given spectra_file file existence and format"""
        folder = self.folder
        data_file = self.data_file + '.pdf'
        spectra_file = self.spectra_file
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_spectra_file_exist(self):
        """assert given spectra_file file existence and format"""
        folder = self.folder
        data_file = self.data_file
        spectra_file = 'file_does_not_exist' + '.txt'
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_data_file_exist(self):
        """assert given spectra_file file existence and format"""
        folder = self.folder
        data_file = 'file_does_not_exist' + '.txt'
        spectra_file = self.spectra_file
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_repeat(self):
        folder = self.folder
        data_file = self.data_file
        spectra_file = self.spectra_file
        repeat = -1
        self.assertRaises(ValueError, Experiment, repeat=repeat, data_file=data_file, spectra_file=spectra_file,
                          folder=folder)
        repeat = 3.6
        self.assertRaises(ValueError, Experiment, repeat=repeat, data_file=data_file, spectra_file=spectra_file,
                          folder=folder)

    def test_load_txt_csv(self):
        experiment = Experiment(data_file='_data_xy_unit_test.txt', spectra_file=self.spectra_file, folder=self.folder)
        _dict_expected = np.array([1.003423, 1.008694, 1.008373, 1.004356, 1.008168, 1.016091])
        _dict_returned = np.array(experiment.data[0])

        self.assertEqual(_dict_returned[0], _dict_expected[0])
        self.assertEqual(_dict_returned[1], _dict_expected[1])
        self.assertEqual(_dict_returned[2], _dict_expected[2])
        self.assertEqual(_dict_returned[3], _dict_expected[3])
        self.assertEqual(_dict_returned[4], _dict_expected[4])
        self.assertEqual(_dict_returned[5], _dict_expected[5])

    def test_loaded_data_sep(self):
        folder = self.folder
        data_file = '_data_sep_unit_test.txt'
        spectra_file = self.spectra_file
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_loaded_data_str(self):
        folder = self.folder
        data_file = '_data_str_unit_test.txt'
        spectra_file = self.spectra_file
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_loaded_spectra_sep(self):
        folder = self.folder
        data_file = self.spectra_file
        spectra_file = '_data_sep_unit_test.txt'
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_loaded_spectra_str(self):
        folder = self.folder
        data_file = self.spectra_file
        spectra_file = '_data_str_unit_test.txt'
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_xy_scaled(self):
        experiment = Experiment(data_file=self.data_file, spectra_file=self.spectra_file, folder=self.folder)
        x_interp, y_interp = experiment.xy_scaled(energy_min=self.energy_min,
                                                  energy_max=self.energy_max,
                                                  energy_step=self.energy_step,
                                                  x_type='energy',
                                                  y_type='attenuation',
                                                  offset_us=0,
                                                  source_to_detector_m=15)

        self.assertAlmostEqual(x_interp[1] - x_interp[0], self.energy_step, delta=self.energy_step / 1000)

    def test_x_raw(self):
        experiment = Experiment(data_file='_data_xy_unit_test.txt', spectra_file=self.spectra_file, folder=self.folder)
        _x_returned = experiment.get_x(x_type='energy', offset_us=0., source_to_detector_m=15)
        _x_expected = np.array([5.825324e+00,
                                5.821177e+00,
                                5.817034e+00,
                                5.812896e+00])
        self.assertAlmostEqual(_x_returned[-1], _x_expected[-1], delta=0.000001)
        self.assertAlmostEqual(_x_returned[-2], _x_expected[-2], delta=0.000001)
        self.assertAlmostEqual(_x_returned[-3], _x_expected[-3], delta=0.000001)
        self.assertAlmostEqual(_x_returned[-4], _x_expected[-4], delta=0.000001)
        _x_returned = experiment.get_x(x_type='lambda', offset_us=0., source_to_detector_m=15)
        _x_expected = np.array([0.118490,
                                0.118532,
                                0.118575,
                                0.118617])
        self.assertAlmostEqual(_x_returned[-1], _x_expected[-1], delta=0.000001)
        self.assertAlmostEqual(_x_returned[-2], _x_expected[-2], delta=0.000001)
        self.assertAlmostEqual(_x_returned[-3], _x_expected[-3], delta=0.000001)
        self.assertAlmostEqual(_x_returned[-4], _x_expected[-4], delta=0.000001)

    def test_y_raw(self):
        experiment = Experiment(data_file='_data_xy_unit_test.txt', spectra_file=self.spectra_file, folder=self.folder)
        _y_returned = experiment.get_y(y_type='transmission', baseline=False)
        _y_expected = np.array([1.003423, 1.008694, 1.008373, 1.004356, 1.008168, 1.016091])
        self.assertAlmostEqual(_y_returned[-1], _y_expected[-1], delta=0.000001)
        self.assertAlmostEqual(_y_returned[-2], _y_expected[-2], delta=0.000001)
        self.assertAlmostEqual(_y_returned[-3], _y_expected[-3], delta=0.000001)
        self.assertAlmostEqual(_y_returned[-4], _y_expected[-4], delta=0.000001)
        _y_returned = experiment.get_y(y_type='attenuation', baseline=False)
        _y_expected = np.array([-0.003423, -0.008694, -0.008373, -0.004356, -0.008168, -0.016091])
        self.assertAlmostEqual(_y_returned[-1], _y_expected[-1], delta=0.000001)
        self.assertAlmostEqual(_y_returned[-2], _y_expected[-2], delta=0.000001)
        self.assertAlmostEqual(_y_returned[-3], _y_expected[-3], delta=0.000001)
        self.assertAlmostEqual(_y_returned[-4], _y_expected[-4], delta=0.000001)
