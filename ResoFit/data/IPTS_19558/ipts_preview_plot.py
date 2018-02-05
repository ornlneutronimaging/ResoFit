from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
import numpy as np
import pprint
import matplotlib.pyplot as plt
from ResoFit.experiment import Experiment
import peakutils as pku

folder = 'data/IPTS_19558/reso_data_19558'
data_file1 = 'Gd_thin.csv'
data_file2 = 'Gd_thick.csv'
data_file3 = 'spheres.csv'
spectra_file = 'Image002_Spectra.txt'

repeat = 1
source_to_detector_m = 16.45  # 16#16.445359069030175#16.447496101100739
offset_us = 2.752  # 0#2.7120797253959119#2.7355447625559037
# transmission = False
baseline = False
energy_xmax = 150
lambda_xmax = None
x_axis = 'number'

# # Calibrate the peak positions
experiment1 = Experiment(data_file=data_file1,
                         spectra_file=spectra_file,
                         repeat=repeat,
                         folder=folder)
experiment2 = Experiment(data_file=data_file2,
                         spectra_file=spectra_file,
                         repeat=repeat,
                         folder=folder)

experiment1.plot_raw(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
                     x_type=x_axis, baseline=baseline, energy_xmax=energy_xmax,
                     lambda_xmax=lambda_xmax)
experiment2.plot_raw(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
                     x_type=x_axis, baseline=baseline, energy_xmax=energy_xmax,
                     lambda_xmax=lambda_xmax)
# experiment1.export_raw(offset_us=offset_us)

# x1, y1 = experiment1.xy_scaled(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
#                                energy_min=7, energy_max=150, energy_step=0.01)
# plt.plot(x1, y1, 'r.', label='interp1')
# x2, y2 = experiment2.xy_scaled(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
#                                energy_min=7, energy_max=150, energy_step=0.01)
# plt.plot(x2, y2, 'k.', label='interp2')
plt.show()
