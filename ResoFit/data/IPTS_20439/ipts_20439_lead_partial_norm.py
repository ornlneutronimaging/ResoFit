import numpy as np
import pprint
import matplotlib.pyplot as plt
from ResoFit.experiment import Experiment
import peakutils as pku
from ResoFit.simulation import Simulation
from scipy import signal
import scipy

folder = 'data/IPTS_20439/reso_data_20439'
sample_name = ['No Pb', '10mm Pb']
data_file = ['Ta_no_lead_partial.csv', 'Ta_10mm_lead_partial.csv']
norm_to_file = ['blank_no_lead_partial.csv', 'blank_no_lead_partial.csv']
norm_factor = [0.95, 0.7]
spectra_file = 'Ta_lead_10mm__0__040_Spectra.txt'


baseline = False
deg = 6
# x_axis = 'number'
logx = False
# # # Calibrate the peak positions
x_type = 'energy'
y_type = 'transmission'
source_to_detector_m = 16.45
offset_us = 0
fmt = '-'
lw = 1

exps = {}
ax0 = None

for _index, each_name in enumerate(sample_name):
    exps[each_name] = Experiment(spectra_file=spectra_file, data_file=data_file[_index], folder=folder)
    exps[each_name].norm_to(file=norm_to_file[_index], norm_factor=norm_factor[_index])
    if ax0 is None:
        ax0 = exps[each_name].plot(x_type=x_type, y_type=y_type,
                                   source_to_detector_m=source_to_detector_m, offset_us=offset_us,
                                   logx=logx, baseline=baseline, deg=deg, fmt=fmt, lw=lw, label=each_name)
    else:
        ax0 = exps[each_name].plot(ax_mpl=ax0, x_type=x_type, y_type=y_type,
                                   source_to_detector_m=source_to_detector_m, offset_us=offset_us,
                                   logx=logx, baseline=baseline, deg=deg, fmt=fmt, lw=lw, label=each_name)

        # simu.plot(ax_mpl=ax0[i], x_type='energy', y_type='attenuation',
        #           source_to_detector_m=source_to_detector_m, offset_us=offset_us, logx=True,
        #           mixed=False, all_layers=False, all_elements=False, items_to_plot=[_ele],
        #           fmt='-.', lw=1, alpha=1)

plt.xlim(5, 120)

plt.show()
#
# # experiment1.plot(offset_us=offset_us, source_to_detector_m=source_to_detector_m,
# #                      x_axis=x_axis, baseline=baseline, energy_xmax=energy_xmax,
# #                      lambda_xmax=lambda_xmax)
#
# data = 1-experiment1.data[0]
# # plt.plot(data, 'k-')
# _baseline = pku.baseline(data, deg=7)
# data_flat = data - _baseline
# plt.plot(data_flat, 'b-')
#
#
# indexes = pku.indexes(data_flat, thres=0.1, min_dist=50)
# print(indexes)
# plt.plot(data_flat[indexes], 'bx', label='peak')
#
# # peakind = signal.find_peaks_cwt(data_flat, widths=np.arange(1, len(data_flat)))
# # print(peakind)
# # plt.plot(data_flat[peakind], 'bs', label='peak')
#
# # After slicing
# experiment1.slice(slice_start=300, slice_end=2200, reset_index=True)
# data_sliced = 1-experiment1.data[0]
# # plt.plot(data_sliced, 'r:')
# _baseline_2 = pku.baseline(data_sliced, deg=7)
# data_sliced_flat = data_sliced - _baseline_2
# plt.plot(experiment1.img_num, data_sliced_flat, 'y-')
#
# indexes = pku.indexes(data_sliced_flat, thres=0.1, min_dist=50)
# x_indexes = indexes + 300
# print(indexes)
# plt.plot(x_indexes, data_sliced_flat[indexes], 'rx', label='peak')
#
# # peakind = signal.find_peaks_cwt(data_sliced_flat, widths=np.arange(1, len(data_sliced_flat)))
# # print(peakind)
# # plt.plot(data_sliced_flat[peakind], 'rs', label='peak')
#
# plt.show()


# plt.plot(x,y)
# plt.show()
# simulation.plot(items_to_plot=['U233'])