from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
import matplotlib.pyplot as plt
import numpy as np
from ResoFit._utilities import get_foil_density_gcm3
from ResoFit._utilities import Layer
import pprint

# Global parameters
energy_min = 4.02
energy_max = 1000
energy_step = 0.01
# Input sample name or names as str, case sensitive
layers = Layer()
# layers.add_layer(layer='Ag', thickness_mm=0.025)
# layers.add_layer(layer='Co', thickness_mm=0.025)
# layers.add_layer(layer='Hf', thickness_mm=0.025)
layers.add_layer(layer='W', thickness_mm=0.05)
# layers.add_layer(layer='In', thickness_mm=0.05)
# layers.add_layer(layer='Cd', thickness_mm=0.5)
# layers.add_layer(layer='Au', thickness_mm=0.01)
# simu = Simulation(energy_min=energy_min, energy_max=energy_max, energy_step=energy_step)
# simu.add_Layer(layer=layers)
# peak_dict = simu.peak_map(thres=0.015, min_dist=20)
# pprint.pprint(peak_dict)


folder = 'data/IPTS_13639/reso_data_13639'
data_file = 'W.csv'
spectra_file = 'spectra.csv'
image_start = 300  # Can be omitted or =None
image_end = 2720  # Can be omitted or =None
# norm_to_file = 'ob_1.csv'  #'Ag.csv'
# norm_to_file = 'Ag.csv'
norm_to_file = 'Hf.csv'
baseline = True
each_step = False

repeat = 0.8
source_to_detector_m = 16.126845685903064  # 16#16.445359069030175#16.447496101100739
offset_us = -12112.431834715671  # 0#2.7120797253959119#2.7355447625559037

# Calibrate the peak positions
calibration = Calibration(data_file=data_file,
                          spectra_file=spectra_file,
                          layer=layers,
                          energy_min=energy_min,
                          energy_max=energy_max,
                          energy_step=energy_step,
                          folder=folder,
                          baseline=baseline)

calibration.experiment.norm_to(norm_to_file, norm_factor=repeat)
calibration.experiment.slice(start=image_start, end=image_end)
# calibration.experiment.plot(
#     x_type='energy',
#     source_to_detector_m=source_to_detector_m,
#     offset_us=offset_us,
#     logx=True,
#     fmt='-'
# )
calibrate_result = calibration.calibrate(source_to_detector_m=source_to_detector_m,
                                         offset_us=offset_us,
                                         vary='none',
                                         each_step=each_step)
calibration.index_peak(thres_exp=0.05, min_dist_exp=5, min_dist_map=5, thres_map=0.05, rel_tol=0.015)
# calibration.analyze_peak()
calibration.plot(y_type='attenuation',
                 # y_type='transmission',
                 x_type='energy',
                 # t_unit='ms',
                 # before=True,
                 # interp=True,
                 mixed=True,
                 # peak_exp='all',
                 table=False,
                 peak_exp='indexed',
                 peak_height=True,
                 index_level='iso',
                 peak_id='all',
                 logx=True,
                 )
plt.xlim(left=0, right=1000)
plt.show()

df = calibration.export(y_type='attenuation',
                        # y_type='transmission',
                        x_type='energy',
                        # t_unit='ms',
                        # before=True,
                        # interp=True,
                        # mixed=True,
                        index_level='iso',
                        peak_id='all')

# # Fit the peak height
# fit = FitResonance(folder=folder,
#                    spectra_file=spectra_file,
#                    data_file=data_file,
#                    repeat=repeat,
#                    energy_min=energy_min,
#                    energy_max=energy_max,
#                    energy_step=energy_step,
#                    calibrated_offset_us=calibration.calibrated_offset_us,
#                    calibrated_source_to_detector_m=calibration.calibrated_source_to_detector_m,
#                    norm_to_file=norm_to_file,
#                    slice_start=image_start,
#                    slice_end=image_end,
#                    baseline=baseline)
# fit_result = fit.fit(layer, vary='thickness', each_step=each_step)
# fit.molar_conc()
# fit.plot()
#
