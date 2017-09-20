import unittest
from ResoFit.experiment import Experiment


class TestExperiment(unittest.TestCase):

    folder = 'data'
    data_file = '_data_unit_test.txt'
    spectra_file = '_spectra_unit_test.txt'

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
        data_file = self.data_file + '.txt'
        spectra_file = self.spectra_file + '.txt'
        repeat = -1
        self.assertRaises(ValueError, Experiment, repeat=repeat, data_file=data_file, spectra_file=spectra_file, folder=folder)
        repeat = 3.5
        self.assertRaises(ValueError, Experiment, repeat=repeat, data_file=data_file, spectra_file=spectra_file, folder=folder)

    # def test_y_raw(self):
    #     pass


        # def test_x(self):
        #     self.fail()
        #
        # def test_ob_y(self):
        #     self.fail()
        #
        # def test_y(self):
        #     self.fail()
