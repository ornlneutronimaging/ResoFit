from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
import numpy as np
import pprint
import matplotlib.pyplot as plt
from ResoFit.experiment import Experiment
import peakutils as pku

folder = 'data'
data_file1 = 'Gd_thin.csv'
data_file2 = 'Gd_thick.csv'
spectra_file = 'Image002_Spectra.txt'

repeat = 1
source_to_detector_m = 16.45  # 16#16.445359069030175#16.447496101100739
offset_us = 2.752  # 0#2.7120797253959119#2.7355447625559037
angstrom = False
transmission = False
baseline = True
energy_xmax = 150
lambda_xmax = None
x_axis = 'time'

# # Calibrate the peak positions
experiment1 = Experiment(data_file=data_file1,
                         spectra_file=spectra_file,
                         repeat=repeat,
                         folder=folder)
experiment2 = Experiment(data_file=data_file2,
                         spectra_file=spectra_file,
                         repeat=repeat,
                         folder=folder)
# exp_x, exp_y = experiment.xy_scaled(energy_min=energy_min,
#                                     energy_max=energy_max,
#                                     energy_step=energy_step,
#                                     offset_us=offset_us,
#                                     source_to_detector_m=source_to_detector_m)
# experiment.norm_to('Ag.csv')
# experiment.slice(slice_start=500, slice_end=1600)

experiment1.plot_raw(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
                     x_axis=x_axis, transmission=transmission, baseline=baseline, energy_xmax=energy_xmax,
                     lambda_xmax=lambda_xmax)
experiment2.plot_raw(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
                     x_axis=x_axis, transmission=transmission, baseline=baseline, energy_xmax=energy_xmax,
                     lambda_xmax=lambda_xmax)
plt.show()
