from unittest import TestCase
from ResoFit.experiment import Experiment

class TestExperiment(TestCase):

    def test_folder(self):
        '''assert given folder existence'''
        folder = 'folder_not_exist'
        data = 'all_thin.txt'
        spectra = 'Image002_Spectra.txt'
        self.assertRaises(ValueError, Experiment, data=data, spectra=spectra, folder=folder)

    def test_spectra_file(self):
        '''assert given spectra file existence and format'''
        folder = 'data'
        data = 'all_thin.txt'
        spectra = 'file_not_exist.txt'
        self.assertRaises(ValueError, Experiment, data=data, spectra=spectra, folder=folder)
        spectra = 'Image002_Spectra.tt'
        self.assertRaises(ValueError, Experiment, data=data, spectra=spectra, folder=folder)

    def test_data_file(self):
        '''assert given data file existence and format'''
        folder = 'data'
        data = 'file_not_exist.txt'
        spectra = 'Image002_Spectra.txt'
        self.assertRaises(ValueError, Experiment, data=data, spectra=spectra, folder=folder)
        data = 'all_thin.tt'
        self.assertRaises(ValueError, Experiment, data=data, spectra=spectra, folder=folder)

    # def test_x(self):
    #     self.fail()
    #
    # def test_ob_y(self):
    #     self.fail()
    #
    # def test_y(self):
    #     self.fail()
