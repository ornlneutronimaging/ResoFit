import unittest
from ResoFit.simulation import Simulation
import numpy as np


class TestCalibration(unittest.TestCase):
    energy_min = 7
    energy_max = 10
    energy_step = 1