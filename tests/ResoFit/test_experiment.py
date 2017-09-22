import unittest
from ResoFit.experiment import Experiment


class TestExperiment(unittest.TestCase):
    folder = 'data/_mock_data_for_test'
    data_file = '_data_unit_test.txt'
    spectra_file = '_spectra_unit_test.txt'
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
        _dict_expected = {0: [1, 2, 3, 4, 5, 6],
                          1: [1.003423, 1.008694, 1.008373, 1.004356, 1.008168, 1.016091]}
        _dict_returned = experiment.data.to_dict('list')

        self.assertEqual(_dict_returned, _dict_expected)

    def test_loaded_data(self):
        folder = self.folder
        data_file = '_data_sep_unit_test.txt'
        spectra_file = self.spectra_file
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_xy_scaled(self):
        experiment = Experiment(data_file=self.data_file, spectra_file=self.spectra_file, folder=self.folder)
        x_interp, y_interp = experiment.xy_scaled(energy_min=self.energy_min,
                                                  energy_max=self.energy_max,
                                                  energy_step=self.energy_step,
                                                  angstrom=False,
                                                  transmission=False,
                                                  offset_us=0,
                                                  source_to_detector_m=15)

        self.assertAlmostEqual(x_interp[1]-x_interp[0], self.energy_step, delta=self.energy_step/1000)


        # def test_y_raw(self):
        #     pass
