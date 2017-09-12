import unittest
import numpy as np
import os
import pprint
from ResoFit._utilities import Experiment


class TestSimulation(unittest.TestCase):




    pass

class TestExperiment(unittest.TestCase):

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

    # def test_file(self):
    #     folder = 'data'
    #     filename = 'Image002_Spectra.txt'
    #     _answer = os.path.exists(folder + '/' + filename)
    #     self.assertTrue(_answer)
    #     filename = 'Image002_Spectr.txt'
    #     _answer = os.path.exists(folder + '/' + filename)
    #     self.assertFalse(_answer)
    #
    # def test_txt(self):
    #     folder = 'data'
    #     filename = 'Image002_Spectra.tt'
    #     self.assertRaises(ValueError, Experiment, filename=filename, folder=folder)