import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation
from ResoFit.calibration import Calibration
import numpy as np
import pprint

# Global parameters
energy_min = 7
energy_max = 100
energy_step = 0.01
# Input sample name or names as str, case sensitive
_layer_1 = 'Gd'
_thickness_1 = 0.15  # mm
mass = 0.36  # gram
length = 25
width = 25
height = 0.025
mm3_to_cm3 = 0.001
density = mass / (length * width * height * mm3_to_cm3)
_density_1 = np.NaN
data_file = 'all_thick.txt'
spectra_file = 'Image002_Spectra.txt'
repeat = 5
source_to_detector_m = 16.45#16#16.445359069030175#16.447496101100739
offset_us = 2.788#0#2.7120797253959119#2.7355447625559037

# # Ideal
# simulation = Simulation(layer_1=_layer_1,
#                         thickness_1=_thickness_1,
#                         density_1=np.NaN,
#                         energy_min=energy_min,
#                         energy_max=energy_max,
#                         energy_step=energy_step)
# simu_x, simu_y = simulation.xy_simu()
#
# plt.plot(simu_x, simu_y, 'b.', label=_layer_1 + '_ideal', markersize=1)
#
# # ideal_y_index = pku.indexes(simu_y, thres=0.15, min_dist=10)  # , thres=0.1, min_dist=50)
# # ideal_x_index = pku.interpolate(simu_x, simu_y, ind=ideal_y_index)
# # print('x_ideal_peak: ', ideal_x_index)
# # plt.plot(simu_x[ideal_y_index], simu_y[ideal_y_index], 'bo', label='peak_ideal')
#
# # Experiment
# experiment = Experiment(data_file=data_file, spectra_file=spectra_file, repeat=repeat)
# exp_x, exp_y = experiment.xy_scaled(energy_min=energy_min,
#                                     energy_max=energy_max,
#                                     energy_step=energy_step,
#                                     source_to_detector_m=source_to_detector_m,
#                                     offset_us=offset_us)
#
# plt.plot(exp_x, exp_y, 'r.', label='data', markersize=1)
#
# # exp_y_index = pku.indexes(exp_y, thres=0.05/max(exp_y), min_dist=7)
# # exp_x_index = pku.interpolate(exp_x, exp_y, ind=exp_y_index)
# # plt.plot(exp_x[exp_y_index], exp_y[exp_y_index], 'go', label='peak_exp')

# Fitting the peak positions
params = Parameters()
params.add('source_to_detector_m', value=source_to_detector_m)
params.add('offset_us', value=offset_us)

calibration = Calibration(data_file=data_file,
                          spectra_file=spectra_file,
                          layer_1=_layer_1,
                          thickness_1=_thickness_1,
                          density_1=np.NaN,
                          energy_min=energy_min,
                          energy_max=energy_max,
                          energy_step=energy_step,
                          repeat=5)

out = calibration.get_exp_params(params_init=params)
calibration.plot()
# out = calibration.exp_params(params)

# out = minimize(peak_y_gap(params), params, method='leastsq',
#                args=(simu_x, simu_y, energy_min, energy_max, energy_step, data, spectra, repeat))
#
# params.pretty_print()
# print(out.__dict__['params'])
# out.__dict__['params'].pretty_print()
# print(out.__dict__['params'].valuesdict()['source_to_detector_m'])
# pprint.pprint(out.__dict__['params'].dumps())

# print(sum((simu_x - exp_x)**2))

# plt.title('Peak estimation')
# plt.ylim(-0.01, 1.01)
# plt.xlim(0, energy_max)
# plt.legend(loc='best')
# plt.show()

#
# df = pd.DataFrame()
# df['Exp_x'] = x_data_array
# df['Exp_y'] = y_data_array
# df2 = pd.DataFrame()
# df2['Ideal_x'] = simu_x
# df2['Ideal_y'] = simu_y
# x_gap = _fit_functions.peak_x_gap(params, ideal_x_index, y_data_array)
# print('x_gap:', x_gap)

# out = minimize(_fit_functions.peak_x_gap, params, method='leastsq', args=(ideal_x_index, y_data_array))
# out = scipy.optimize.minimize(_fit_funtions.peak_x_gap_scipy, delay_us, method='leastsq', args=(ideal_x_index, y_data_array))
# print(out.__dict__)
