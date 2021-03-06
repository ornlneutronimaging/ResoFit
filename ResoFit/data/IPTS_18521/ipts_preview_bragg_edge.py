from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
import numpy as np
import pprint
import matplotlib.pyplot as plt
from ResoFit.experiment import Experiment
import peakutils as pku

folder = 'data/IPTS_18521/bragg_data_18521'
data_file1 = 'run_42_image0_small.txt'
ob_file1 = 'run_42_image0_small_ob.txt'
spectra_file = 'Image017_Spectra.txt'
data_file2 = 'run_42_image1_big_bottom1.txt'
ob_file2 = 'run_42_image1_big_bottom1_ob.txt'
data_file3 = 'run_42_image1_big_bottom2.txt'
ob_file3 = 'run_42_image1_big_bottom2_ob_lines.txt'
data_file4 = 'run_42_image0_big_bottom.txt'
ob_file4 = 'run_42_image0_big_bottom_ob.txt'
repeat = 1
source_to_detector_m = 15.  # 16#16.445359069030175#16.447496101100739
offset_us = 0  # 0#2.7120797253959119#2.7355447625559037
# transmission = False
baseline = False
energy_xmax = 150
lambda_xmax = None
x_axis = 'lambda'

# # Calibrate the peak positions
experiment1 = Experiment(data_file=data_file1,
                         spectra_file=spectra_file,
                         norm_factor=repeat,
                         folder=folder)
experiment1.norm_to(ob_file1)
experiment1.plot(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
                 x_axis=x_axis, baseline=baseline, energy_xmax=energy_xmax,
                 lambda_xmax=lambda_xmax, transmission=True)

experiment2 = Experiment(data_file=data_file2,
                         spectra_file=spectra_file,
                         norm_factor=repeat,
                         folder=folder)
experiment2.norm_to(ob_file2)
experiment2.plot(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
                 x_axis=x_axis, baseline=baseline, energy_xmax=energy_xmax,
                 lambda_xmax=lambda_xmax, transmission=True)

experiment3 = Experiment(data_file=data_file3,
                         spectra_file=spectra_file,
                         norm_factor=repeat,
                         folder=folder)
experiment3.norm_to(ob_file3)
experiment3.plot(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
                 x_axis=x_axis, baseline=baseline, energy_xmax=energy_xmax,
                 lambda_xmax=lambda_xmax, transmission=True)

experiment4 = Experiment(data_file=data_file4,
                         spectra_file=spectra_file,
                         norm_factor=repeat,
                         folder=folder)
experiment4.norm_to(ob_file4)
experiment4.plot(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
                 x_axis=x_axis, baseline=baseline, energy_xmax=energy_xmax,
                 lambda_xmax=lambda_xmax, transmission=True)

# x1, y1 = experiment1.xy_scaled(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
#                                energy_min=7, energy_max=150, energy_step=0.01)
# plt.plot(x1, y1, 'r.', label='interp1')
# x2, y2 = experiment2.xy_scaled(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
#                                energy_min=7, energy_max=150, energy_step=0.01)
# plt.plot(x2, y2, 'k.', label='interp2')
plt.title('IPTS-18521 CT Bragg-edge data of spheres')
plt.show()
