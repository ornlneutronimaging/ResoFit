import unittest
import numpy as np
import os
import pprint
from ResoFit._utilities import Experiment

class TestSpectra(unittest.TestCase):

    def test_filename(self):
        '''assert given folder existence'''
        _file_path = os.path.abspath(os.path.dirname(__file__))
        folder = '/data'
        path_full = _file_path+folder
        _answer = os.path.isdir(path_full)
        self.assertTrue(_answer)
        folder = '/xyz'
        path_full = _file_path+folder
        _answer = os.path.isdir(path_full)
        self.assertFalse(_answer)

    def test_file(self):
        folder = 'data'
        filename = 'Image002_Spectra.txt'
        _answer = os.path.exists(folder + '/' + filename)
        self.assertTrue(_answer)
        filename = 'Image002_Spectr.txt'
        _answer = os.path.exists(folder + '/' + filename)
        self.assertFalse(_answer)

    def test_txt(self):
        folder = 'data'
        filename = 'Image002_Spectra.tt'
        self.assertRaises(ValueError, Experiment, filename=filename, folder=folder)