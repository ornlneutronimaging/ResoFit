import unittest
import numpy as np
import os
import pytest
import ResoFit._utilities as fit_util
from ResoFit.simulation import Simulation


class TestLayer(unittest.TestCase):
    def test_input_type(self):
        layer = fit_util.Layer()
        thickness_1 = 0.018
        density_1 = None
        layer_1 = 16
        pytest.raises(ValueError, layer.add_layer, layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
        # self.assertRaises(ValueError, layer.add_Layer, layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
        layer_1 = []
        pytest.raises(ValueError, layer.add_layer, layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
        layer_1 = {}
        pytest.raises(ValueError, layer.add_layer, layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
        layer_1 = np.array
        pytest.raises(ValueError, layer.add_layer, layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
        layer_1 = 'UO3'
        pytest.raises(ValueError, layer.add_layer, layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
        layer_1 = 'U'
        thickness_1 = ''
        pytest.raises(ValueError, layer.add_layer, layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)

    def test_add_layer(self):
        layer_1 = 'Gd'
        thickness_1 = 0.018
        density_1 = None
        layer = fit_util.Layer()
        layer.add_layer(layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
        info_expected = {'Gd': {'density': {'units': 'g/cm3',
                                            'value': np.NaN,
                                            },
                                'layer': 'Gd',
                                'thickness': {'units': 'mm',
                                              'value': 0.018,
                                              },
                                'molar_mass': {'value': None,
                                               'units': None,
                                               },
                                'molar_conc': {'value': None,
                                               'units': None,
                                               },
                                }
                         }
        assert layer.info == info_expected
        density_1 = 7.9
        layer = fit_util.Layer()
        layer.add_layer(layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
        info_expected = {'Gd': {'density': {'units': 'g/cm3',
                                            'value': 7.9,
                                            },
                                'layer': 'Gd',
                                'thickness': {'units': 'mm',
                                              'value': 0.018,
                                              },
                                'molar_mass': {'value': None,
                                               'units': None,
                                               },
                                'molar_conc': {'value': None,
                                               'units': None,
                                               },
                                }
                         }
        assert layer.info == info_expected


class TestItems(unittest.TestCase):
    layer_1 = 'U'
    thickness_1 = 0.05
    layer_2 = 'Ag'
    thickness_2 = 0.05
    layers = fit_util.Layer()
    layers.add_layer(layer=layer_1, thickness_mm=thickness_1)
    layers.add_layer(layer=layer_2, thickness_mm=thickness_2)
    database = '_data_for_unittest'
    simulation = Simulation(database=database)
    simulation.add_Layer(layers)
    items = fit_util.Items(simulation.o_reso)

    def test_raises(self):
        name = 'AG'
        pytest.raises(ValueError, fit_util._shape_items, name=name)
        name = 'aG'
        pytest.raises(ValueError, fit_util._shape_items, name=name)
        name = 'AgO'
        pytest.raises(ValueError, fit_util._shape_items, name=name)
        name = 'ag'
        pytest.raises(ValueError, fit_util._shape_items, name=name)
        name = ''
        pytest.raises(ValueError, fit_util._shape_items, name=name)
        name = []
        pytest.raises(ValueError, fit_util._shape_items, name=name)

    def test_isotope_format(self):
        name = '238-U'
        expected_path = ['U', 'U', '238-U']
        assert fit_util._shape_items(name) == expected_path
        name = '238U'
        assert fit_util._shape_items(name) == expected_path
        name = 'U-238'
        assert fit_util._shape_items(name) == expected_path
        name = 'U238'
        assert fit_util._shape_items(name) == expected_path

    def test_fill_iso_to_items(self):
        name = 'U*'
        expected_path_list = [['U', 'U', '233-U'],
                              ['U', 'U', '234-U'],
                              ['U', 'U', '235-U'],
                              ['U', 'U', '238-U']]
        assert fit_util._fill_iso_to_items(name, database=self.database) == expected_path_list
        name = 'U'
        pytest.raises(ValueError, fit_util._fill_iso_to_items, name=name, database=self.database)

    def test_shape_items(self):
        name = 'U'
        expected_path = ['U', 'U']
        assert fit_util._shape_items(name) == expected_path
        name = 'u'
        expected_path = ['U', 'U']
        assert fit_util._shape_items(name) == expected_path
        name = 'Gd'
        expected_path = ['Gd', 'Gd']
        assert fit_util._shape_items(name) == expected_path

    def test_items_shaped(self):
        _input = ['Gd', ['U'], 'U-238', 'U*']
        expected = [['Gd', 'Gd'],
                    ['U', 'U'],
                    ['U', 'U', '233-U'],
                    ['U', 'U', '234-U'],
                    ['U', 'U', '235-U'],
                    ['U', 'U', '238-U']]
        obtained = self.items.shaped(_input)
        assert obtained == expected

    def test_items_original(self):
        _input = [['Gd', 'Gd'],
                  ['U', 'U'],
                  ['U', 'U', '233-U'],
                  ['U', 'U', '234-U'],
                  ['U', 'U', '235-U'],
                  ['U', 'U', '238-U']]
        expected = [['Gd', 'Gd'],
                    ['U', 'U'],
                    ['U', 'U', '233-U'],
                    ['U', 'U', '234-U'],
                    ['U', 'U', '235-U'],
                    ['U', 'U', '238-U']]
        obtained = self.items.shaped(_input)
        assert obtained == expected


class TestPeaks(unittest.TestCase):
    energy_min = 7
    energy_max = 150
    energy_step = 0.01
    simulation = Simulation(energy_min=energy_min,
                            energy_max=energy_max,
                            energy_step=energy_step,
                            database='_data_for_unittest')
    simulation.add_layer(layer='U', thickness_mm=0.05)
    x = simulation.o_reso.stack_sigma['U']['U']['energy_eV']
    y = simulation.o_reso.stack_sigma['U']['U']['sigma_b']

    def test_findpeak1(self):
        peak_found = fit_util._find_peak(y=self.y, x=self.x, thres=0.015, min_dist=1, imprv_reso=False)
        peak_expected = {
            'x': [20.87, 36.68, 66.03, 80.75, 102.57, 116.91],
            'y': [9801.184720322835, 13337.61249582717, 4356.4307835150385, 276.22478464487637,
                  6022.958717161858, 2003.9245670422251],
        }
        self.assertDictEqual(peak_found, peak_expected)

    def test_findpeak2(self):
        peak_found = fit_util._find_peak(y=self.y, x=self.x, thres=0.015, min_dist=1, imprv_reso=True)
        peak_expected = {
            'x': [20.87274280816388, 36.68474982964769, 66.03358430830164, 80.75548785869171, 102.56856740218703,
                  116.91012795048718],
            'y': [9801.184720322835, 13337.61249582717, 4356.4307835150385, 276.22478464487637, 6022.958717161858,
                  2003.9245670422251]
        }
        self.assertDictEqual(peak_found, peak_expected)

    def test_findpeak3(self):
        peak_found = fit_util._find_peak(y=self.y, thres=0.015, min_dist=1, imprv_reso=False)
        peak_expected = {
            'x': [1387, 2968, 5903, 7375, 9557, 10991],
            'y': [9801.184720322835, 13337.61249582717, 4356.4307835150385, 276.22478464487637, 6022.958717161858,
                  2003.9245670422251]
        }
        self.assertDictEqual(peak_found, peak_expected)

    def test_findpeak4(self):
        peak_found = fit_util._find_peak(y=self.y, thres=0.015, min_dist=1, imprv_reso=True)
        peak_expected = {
            'x': [1387.2742808163878, 2968.47498296478, 5903.358430830144, 7375.5487858692095, 9556.85674021868,
                  10991.012795048715],
            'y': [9801.184720322835, 13337.61249582717, 4356.4307835150385, 276.22478464487637, 6022.958717161858,
                  2003.9245670422251]
        }
        print(peak_found)
        self.assertDictEqual(peak_found, peak_expected)

    def test_indexes(self):
        x = self.simulation.o_reso.stack_sigma['U']['U']['energy_eV']
        y = self.simulation.o_reso.stack_sigma['U']['U']['sigma_b']
        peak_df = fit_util.find_peak(y=y, x=x, x_name='x', thres=0.015, min_dist=1)
        peak_df_expected = {'x': [20.87, 36.68, 66.03, 80.75, 102.57, 116.91],
                            'y': [9801.18472032, 13337.61249583, 4356.43078352, 276.22478464,
                                  6022.95871716, 2003.92456704],
                            }
        assert peak_df['x'].tolist() == pytest.approx(peak_df_expected['x'])
        assert peak_df['y'].tolist() == pytest.approx(peak_df_expected['y'])

        peak_df = fit_util.find_peak(y=y, x_name='x', thres=0.015, min_dist=1)
        peak_df_expected = {'x': [1387, 2968, 5903, 7375, 9557, 10991],
                            'y': [9801.18472032, 13337.61249583, 4356.43078352, 276.22478464,
                                  6022.95871716, 2003.92456704],
                            }
        assert peak_df['x'].tolist() == pytest.approx(peak_df_expected['x'])
        assert peak_df['y'].tolist() == pytest.approx(peak_df_expected['y'])
        # assert peak_dict == pytest.approx(peak_dict_expected)


def test_get_foil_density_gcm3(length_mm=25, width_mm=25, thickness_mm=0.025, mass_g=0.14):
    expected = 8.96
    assert fit_util.get_foil_density_gcm3(length_mm, width_mm, thickness_mm, mass_g) == expected
