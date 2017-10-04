import unittest
from ResoFit.simulation import Simulation
import numpy as np


class TestSimulation(unittest.TestCase):
    energy_min = 7
    energy_max = 10
    energy_step = 1
    simulation = Simulation(energy_min=energy_min,
                            energy_max=energy_max,
                            energy_step=energy_step)

    def test_add_layer(self):
        simulation = self.simulation
        self.assertIsNone(simulation.simu_x)
        self.assertIsNone(simulation.simu_y)
        simulation.add_layer(layer='U', layer_thickness_mm=0.15)
        _simu_x_returned = simulation.simu_x
        _simu_y_returned = simulation.simu_y
        _simu_x_expected = np.array([7., 8., 9., 10.])
        _simu_y_expected = np.array([0.03699373, 0.00936537, 0.00854215, 0.00726004])
        self.assertAlmostEqual(_simu_x_returned[0], _simu_x_expected[0], delta=0.000001)
        self.assertAlmostEqual(_simu_x_returned[1], _simu_x_expected[1], delta=0.000001)
        self.assertAlmostEqual(_simu_x_returned[2], _simu_x_expected[2], delta=0.000001)
        self.assertAlmostEqual(_simu_x_returned[3], _simu_x_expected[3], delta=0.000001)
        self.assertAlmostEqual(_simu_y_returned[0], _simu_y_expected[0], delta=0.000001)
        self.assertAlmostEqual(_simu_y_returned[1], _simu_y_expected[1], delta=0.000001)
        self.assertAlmostEqual(_simu_y_returned[2], _simu_y_expected[2], delta=0.000001)
        self.assertAlmostEqual(_simu_y_returned[3], _simu_y_expected[3], delta=0.000001)

    def test_set_isotopic_ratio(self):
        simulation = self.simulation
        simulation.add_layer(layer='U', layer_thickness_mm=0.15)
        simulation.set_isotopic_ratio('U', 'U', [0., 0., 0.99, 0.01])
        _isotopic_ratio_list_wrong_len = [0., 0.99, 0.01]
        self.assertRaises(ValueError, simulation.set_isotopic_ratio,
                          layer='U', element='U', new_isotopic_ratio_list=_isotopic_ratio_list_wrong_len)
        _simu_x_returned = simulation.simu_x
        _simu_y_returned = simulation.simu_y
        _simu_x_expected = np.array([7., 8., 9., 10.])
        _simu_y_expected = np.array([0.06464851, 0.01259978, 0.11890677, 0.02255858])
        self.assertAlmostEqual(_simu_x_returned[0], _simu_x_expected[0], delta=0.000001)
        self.assertAlmostEqual(_simu_x_returned[1], _simu_x_expected[1], delta=0.000001)
        self.assertAlmostEqual(_simu_x_returned[2], _simu_x_expected[2], delta=0.000001)
        self.assertAlmostEqual(_simu_x_returned[3], _simu_x_expected[3], delta=0.000001)
        self.assertAlmostEqual(_simu_y_returned[0], _simu_y_expected[0], delta=0.000001)
        self.assertAlmostEqual(_simu_y_returned[1], _simu_y_expected[1], delta=0.000001)
        self.assertAlmostEqual(_simu_y_returned[2], _simu_y_expected[2], delta=0.000001)
        self.assertAlmostEqual(_simu_y_returned[3], _simu_y_expected[3], delta=0.000001)

    def test_x_angstrom(self):
        simulation = self.simulation
        simulation.add_layer(layer='U', layer_thickness_mm=0.15)
        _x_returned = simulation.x_angstrom()
        _x_expected = np.array([0.10809189, 0.10111071, 0.09532809, 0.09043617])
        self.assertAlmostEqual(_x_returned[0], _x_expected[0], delta=0.000001)
        self.assertAlmostEqual(_x_returned[1], _x_expected[1], delta=0.000001)
        self.assertAlmostEqual(_x_returned[2], _x_expected[2], delta=0.000001)
        self.assertAlmostEqual(_x_returned[3], _x_expected[3], delta=0.000001)

    def test_y_tansmission(self):
        simulation = self.simulation
        simulation.add_layer(layer='U', layer_thickness_mm=0.15)
        _y_returned = simulation.y_transmission()
        _y_expected = np.array([ 0.96300627,  0.99063463,  0.99145785,  0.99273996])
        self.assertAlmostEqual(_y_returned[0], _y_expected[0], delta=0.000001)
        self.assertAlmostEqual(_y_returned[1], _y_expected[1], delta=0.000001)
        self.assertAlmostEqual(_y_returned[2], _y_expected[2], delta=0.000001)
        self.assertAlmostEqual(_y_returned[3], _y_expected[3], delta=0.000001)

    def test_xy_simu(self):
        simulation = self.simulation
        simulation.add_layer(layer='U', layer_thickness_mm=0.15)
        _x_returned, _y_returned = simulation.xy_simu(angstrom=True, transmission=True)
        _x_expected = np.array([0.10809189, 0.10111071, 0.09532809, 0.09043617])
        self.assertAlmostEqual(_x_returned[0], _x_expected[0], delta=0.000001)
        self.assertAlmostEqual(_x_returned[1], _x_expected[1], delta=0.000001)
        self.assertAlmostEqual(_x_returned[2], _x_expected[2], delta=0.000001)
        self.assertAlmostEqual(_x_returned[3], _x_expected[3], delta=0.000001)
        _y_expected = np.array([0.96300627, 0.99063463, 0.99145785, 0.99273996])
        self.assertAlmostEqual(_y_returned[0], _y_expected[0], delta=0.000001)
        self.assertAlmostEqual(_y_returned[1], _y_expected[1], delta=0.000001)
        self.assertAlmostEqual(_y_returned[2], _y_expected[2], delta=0.000001)
        self.assertAlmostEqual(_y_returned[3], _y_expected[3], delta=0.000001)
        _x_returned, _y_returned = simulation.xy_simu(angstrom=False, transmission=False)
        _x_expected = np.array([7., 8., 9., 10.])
        self.assertAlmostEqual(_x_returned[0], _x_expected[0], delta=0.000001)
        self.assertAlmostEqual(_x_returned[1], _x_expected[1], delta=0.000001)
        self.assertAlmostEqual(_x_returned[2], _x_expected[2], delta=0.000001)
        self.assertAlmostEqual(_x_returned[3], _x_expected[3], delta=0.000001)
        _y_expected = np.array([0.03699373, 0.00936537, 0.00854215, 0.00726004])
        self.assertAlmostEqual(_y_returned[0], _y_expected[0], delta=0.000001)
        self.assertAlmostEqual(_y_returned[1], _y_expected[1], delta=0.000001)
        self.assertAlmostEqual(_y_returned[2], _y_expected[2], delta=0.000001)
        self.assertAlmostEqual(_y_returned[3], _y_expected[3], delta=0.000001)





