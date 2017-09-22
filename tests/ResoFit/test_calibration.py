import unittest
from ResoFit.calibration import Calibration
import numpy as np


class TestCalibration(unittest.TestCase):
    energy_min = 7
    energy_max = 10
    energy_step = 1

    # def test_vary_tag(self):
    #     calibration = Calibration
