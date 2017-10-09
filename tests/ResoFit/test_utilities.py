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