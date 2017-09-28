from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
from ResoFit.experiment import Experiment
import matplotlib.pyplot as plt
import numpy as np
import pprint

# Global parameters
energy_min = 7
energy_max = 300
energy_step = 0.01
# Input sample name or names as str, case sensitive
layer_1 = 'Co'
thickness_1 = 0.15  # mm
mass = 0.36  # gram
length = 25
width = 25
height = 0.075
mm3_to_cm3 = 0.001
density = np.NaN  # mass / (length * width * height * mm3_to_cm3)

folder = 'data'
data_file = 'Co.csv'
spectra_file = 'spectra.csv'

repeat = 1
source_to_detector_m = 16.123278721983177  # 16#16.445359069030175#16.447496101100739
offset_us = -12112.494119089204  # 0#2.7120797253959119#2.7355447625559037

# # Calibrate the peak positions
# experiment = Experiment(data_file=data_file,
#                         spectra_file=spectra_file,
#                         repeat=repeat,
#                         folder=folder)
# # exp_x, exp_y = experiment.xy_scaled(energy_min=energy_min,
# #                                     energy_max=energy_max,
# #                                     energy_step=energy_step,
# #                                     offset_us=offset_us,
# #                                     source_to_detector_m=source_to_detector_m)
# experiment.norm_to('Ag.csv', baseline=False)
# exp_x_sliced, exp_y_sliced = experiment.slice(400, 2700)
# plt.plot(exp_y_sliced)
# plt.show()

# Calibrate the peak positions
calibration = Calibration(data_file=data_file,
                          spectra_file=spectra_file,
                          layer_1=layer_1,
                          thickness_1=thickness_1,
                          density_1=np.NaN,
                          energy_min=energy_min,
                          energy_max=energy_max,
                          energy_step=energy_step,
                          repeat=repeat,
                          folder=folder)

calibration.norm_to('Ag.csv')
# calibration.slice(400, 2400)


calibrate_result = calibration.calibrate(source_to_detector_m=source_to_detector_m,
                                         offset_us=offset_us,
                                         vary='offset')
calibration.plot_before()
calibration.plot_after()
# calibration.plot_after_interp()

# Fit the peak height
fit = FitResonance(spectra_file=spectra_file,
                   data_file=data_file,
                   repeat=repeat,
                   layer=layer_1,
                   energy_min=energy_min,
                   energy_max=energy_max,
                   energy_step=energy_step,
                   calibrated_offset_us=calibration.calibrated_offset_us,
                   calibrated_source_to_detector_m=calibration.calibrated_source_to_detector_m,
                   norm_to_file='Ag.csv')#,
                   # slice_start=400,
                   # slice_end=2400)
fit.fit(thickness=thickness_1, density=density, vary='thickness')
fit.molar_conc(layer_1)
fit.plot_before()
fit.plot_after()
