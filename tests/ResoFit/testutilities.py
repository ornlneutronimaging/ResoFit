from unittest import TestCase
import pandas as pd


class TestLoad(TestCase):
    def test_load_txt_csv(self):
        from ResoFit._utilities import load_txt_csv
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        y = [1.003423,
             1.008694,
             1.008373,
             1.004356,
             1.008168,
             1.016091,
             1.019222,
             1.015934,
             1.012372]
        df_expected = pd.DataFrame(y, index=None)
        df_expected.insert(0, x)
        df = load_txt_csv('ResoFit/data/_data_xy_unit_test.txt')

        self.assertEqual(df, df_expected)
