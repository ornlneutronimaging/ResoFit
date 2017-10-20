import unittest
import numpy as np
import os
import ResoFit._utilities as fit_util


class TestLayer(unittest.TestCase):
    def test_layer_type(self):
        layer_1 = 16
        thickness_1 = 0.018
        density_1 = None
        layer = fit_util.Layer()
        self.assertRaises(ValueError, layer.add_layer, layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
        layer_1 = []
        self.assertRaises(ValueError, layer.add_layer, layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
        layer_1 = {}
        self.assertRaises(ValueError, layer.add_layer, layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)

    def test_add_layer(self):
        layer_1 = 'Gd'
        thickness_1 = 0.018
        density_1 = None
        layer = fit_util.Layer()
        layer.add_layer(layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
        info_expected = {'Gd': {'density': {'units': 'g/cm3',
                                            'value': np.NaN},
                                'layer': 'Gd',
                                'thickness': {'units': 'mm',
                                              'value': 0.018},
                                'molar_mass': {'value': None,
                                               'units': None},
                                'molar_conc': {'value': None,
                                               'units': None}
                                }
                         }
        self.assertEqual(layer.info, info_expected)
        density_1 = 7.9
        layer = fit_util.Layer()
        layer.add_layer(layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
        info_expected = {'Gd': {'density': {'units': 'g/cm3',
                                            'value': 7.9},
                                'layer': 'Gd',
                                'thickness': {'units': 'mm',
                                              'value': 0.018},
                                'molar_mass': {'value': None,
                                               'units': None},
                                'molar_conc': {'value': None,
                                               'units': None}
                                }
                         }
        self.assertEqual(layer.info, info_expected)


class TestRestructureInput(unittest.TestCase):
    def test_element(self):
        expected_path = ['Gd', 'Gd']
        name = 'Gd'
        assert fit_util.restructure_input(name) == expected_path
        name = 'GD'
        self.assertRaises(ValueError, fit_util.restructure_input, name=name)
        name = 'gD'
        self.assertRaises(ValueError, fit_util.restructure_input, name=name)
        name = 'GdO'
        self.assertRaises(ValueError, fit_util.restructure_input, name=name)
        name = 'gd'
        self.assertRaises(ValueError, fit_util.restructure_input, name=name)

    def test_isotope_format_1(self):
        name = '238-U'
        expected_path = ['U', 'U', '238-U']
        assert fit_util.restructure_input(name) == expected_path

    # def test_isotope_format_2(self):
    #     name = '238U'
    #     expected_path = ['U', 'U', '238-U']
    #     assert fit_util.restructure_input(name) == expected_path
    #
    # def test_isotope_format_3(self):
    #     name = 'U-238'
    #     expected_path = ['U', 'U', '238-U']
    #     assert fit_util.restructure_input(name) == expected_path
    #
    # def test_isotope_format_4(self):
    #     name = 'U238'
    #     expected_path = ['U', 'U', '238-U']
    #     assert fit_util.restructure_input(name) == expected_path