from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
import matplotlib.pyplot as plt
import numpy as np
from ResoFit._utilities import get_foil_density_gcm3
from ResoFit._utilities import Layer
import pprint

# Global parameters
energy_min = 4.09
energy_max = 1000
energy_step = 0.01
# Input sample name or names as str, case sensitive
layers = Layer()
layers.add_layer(layer='Ag', thickness_mm=0.025)
# layers.add_layer(layer='Co', thickness_mm=0.025)
# layers.add_layer(layer='Hf', thickness_mm=0.025)
# layers.add_layer(layer='W', thickness_mm=0.05)
# layers.add_layer(layer='In', thickness_mm=0.05)
# layers.add_layer(layer='Cd', thickness_mm=0.5)
# layers.add_layer(layer='Au', thickness_mm=0.01)
# simu = Simulation(energy_min=energy_min, energy_max=energy_max, energy_step=energy_step)
# simu.add_Layer(layer=layers)
# peak_dict = simu.peak_map(thres=0.015, min_dist=20)
# pprint.pprint(peak_dict)


folder = 'data/IPTS_13639/reso_data_13639'
data_file = 'Ag.csv'
spectra_file = 'spectra.csv'
image_start = 300  # Can be omitted or =None
image_end = 2720  # Can be omitted or =None
# norm_to_file = 'ob_1.csv'  #'Ag.csv'
# norm_to_file = 'Ag.csv'
norm_to_file = 'ob_all.csv'
baseline = False
each_step = True

repeat = 1.2
source_to_detector_m = 16.126845685903064  # 16#16.445359069030175#16.447496101100739
offset_us = -12112.431834715671  # 0#2.7120797253959119#2.7355447625559037

# Calibrate the peak positions
calibration = Calibration(data_file=data_file,
                          spectra_file=spectra_file,
                          layer=layers,
                          energy_min=energy_min,
                          energy_max=energy_max,
                          energy_step=energy_step,
                          norm_factor=repeat,
                          folder=folder,
                          baseline=baseline)

calibration.experiment.norm_to(norm_to_file)
calibration.experiment.slice(start=image_start, end=image_end)

calibration.experiment.plot(
    x_type='energy',
    source_to_detector_m=source_to_detector_m,
    offset_us=offset_us,
    logx=True,
)
# calibrate_result = calibration.calibrate(source_to_detector_m=source_to_detector_m,
#                                          offset_us=offset_us,
#                                          vary='all',
#                                          each_step=each_step)
# calibration.index_peak(thres=0.05, min_dist=10, map_min_dist=10, map_thres=0.05)
# # calibration.analyze_peak()
# calibration.plot(
#     y_type='attenuation',
#     # y_type='transmission',
#     x_type='energy',
#     # t_unit='ms',
#     # before=True,
#     # interp=True,
#     mixed=True,
#     # peak_exp='all',
#     table=False,
#     peak_exp='indexed',
#     peak_height=True,
#     index_level='iso',
#     peak_id='all',
#     logx=True,
# )
plt.xlim(left=0, right=1000)
plt.show()
