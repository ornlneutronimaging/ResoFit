import numpy as np
import pprint
import matplotlib.pyplot as plt
from ResoFit.experiment import Experiment
import peakutils as pku
# from scipy import signal

folder = 'data/IPTS_19558/reso_data_19558'
data_file1 = 'spheres.csv'
# data_file2 = 'spheres_background_1.csv'
spectra_file = 'Image002_Spectra.txt'

source_to_detector_m = 16.45  # 16#16.445359069030175#16.447496101100739
offset_us = 2.752  # 0#2.7120797253959119#2.7355447625559037
baseline = False
energy_xmax = 150
lambda_xmax = None
x_axis = 'number'

# # Calibrate the peak positions
experiment1 = Experiment(data_file=data_file1,
                         spectra_file=spectra_file,
                         folder=folder)

# experiment1.plot_raw(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
#                      x_axis=x_axis, baseline=baseline, energy_xmax=energy_xmax,
#                      lambda_xmax=lambda_xmax)

data = 1-experiment1.data[0]
# plt.plot(data, 'k-')
_baseline = pku.baseline(data, deg=7)
data_flat = data - _baseline
plt.plot(data_flat, 'b-')


indexes = pku.indexes(data_flat, thres=0.1, min_dist=50)
print(indexes)
plt.plot(data_flat[indexes], 'bx', label='peak')

# peakind = signal.find_peaks_cwt(data_flat, widths=np.arange(1, len(data_flat)))
# print(peakind)
# plt.plot(data_flat[peakind], 'bs', label='peak')

# After slicing
experiment1.slice(slice_start=300, slice_end=2200, reset_index=True)
data_sliced = 1-experiment1.data[0]
# plt.plot(data_sliced, 'r:')
_baseline_2 = pku.baseline(data_sliced, deg=7)
data_sliced_flat = data_sliced - _baseline_2
plt.plot(experiment1.img_num, data_sliced_flat, 'y-')

indexes = pku.indexes(data_sliced_flat, thres=0.1, min_dist=50)
x_indexes = indexes + 300
print(indexes)
plt.plot(x_indexes, data_sliced_flat[indexes], 'rx', label='peak')

# peakind = signal.find_peaks_cwt(data_sliced_flat, widths=np.arange(1, len(data_sliced_flat)))
# print(peakind)
# plt.plot(data_sliced_flat[peakind], 'rs', label='peak')

plt.show()
# experiment2 = Experiment(data_file=data_file2,
#                          spectra_file=spectra_file,
#                          repeat=repeat,
#                          folder=folder)
# experiment2.plot_raw(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
#                      x_axis=x_axis, baseline=baseline, energy_xmax=energy_xmax,
#                      lambda_xmax=lambda_xmax)


# experiment1.export_raw(offset_us=offset_us)


plt.show()
