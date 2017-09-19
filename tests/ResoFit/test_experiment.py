import unittest
from ResoFit.experiment import Experiment


class TestExperiment(unittest.TestCase):
    def test_folder(self):
        """assert given folder existence"""
        folder = 'data'
        data_file = 'all_thin.txt'
        spectra_file = 'Image002_spectra_file.txt'
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_spectra_file_exist(self):
        """assert given spectra_file file existence and format"""
        folder = 'data'
        data_file = 'all_thin.txt'
        spectra_file = 'Image002_spectra_file.txt'

        # spectra_file = 'file_does_not_exist.txt'
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_spectra_file_format(self):
        """assert given spectra_file file existence and format"""
        folder = 'data'
        data_file = 'all_thin.txt'
        spectra_file = 'Image002_spectra_file.txt'
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)
        spectra_file = 'Image002_spectra_file.tt'
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_data_file(self):
        """assert given data_file file existence and format"""
        folder = 'data'
        data_file = 'file_not_exist.txt'
        spectra_file = 'Image002_spectra_file.txt'
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)
        data_file = 'all_thin.tt'
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)

    def test_repeat(self):
        folder = 'data'
        data_file = 'all_thin.txt'
        spectra_file = 'file_not_exist.txt'
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)
        spectra_file = 'Image002_spectra_file.tt'
        self.assertRaises(ValueError, Experiment, data_file=data_file, spectra_file=spectra_file, folder=folder)


        # def test_x(self):
        #     self.fail()
        #
        # def test_ob_y(self):
        #     self.fail()
        #
        # def test_y(self):
        #     self.fail()
