from unittest import TestCase
import pandas as pd


class TestLoad_txt_csv(TestCase):
    def test_load_txt_csv(self):
        from ResoFit._utilities import load_txt_csv
        _dict_expected = {0: [1, 2, 3, 4, 5, 6],
                          1: [1.003423, 1.008694, 1.008373, 1.004356, 1.008168, 1.016091]}
        # df_expected.insert(1, column=1, value=y)
        _df = load_txt_csv('_data_xy_unit_test.txt')
        _dict_returned = _df.to_dict('list')
        print(_dict_returned)
        self.assertEqual(_dict_returned, _dict_expected)

