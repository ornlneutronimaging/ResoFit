from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation
import matplotlib.pyplot as plt
import numpy as np
import pprint
from ResoFit._utilities import get_foil_density_gcm3
from ResoFit._utilities import Layer
import pprint


# Global parameters
energy_min = 14
energy_max = 300
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
image_start = 500  # Can be omitted or =None
image_end = 1600  # Can be omitted or =None
# norm_to_file = 'ob_1.csv'  #'Ag.csv'
# norm_to_file = 'Ag.csv'
norm_to_file = 'ob_all.csv'
baseline = True
each_step = False

repeat = 1
source_to_detector_m = 16.123278721983177  # 16#16.445359069030175#16.447496101100739
offset_us = -12112.494119089204  # 0#2.7120797253959119#2.7355447625559037

# Calibrate the peak positions
calibration = Calibration(data_file=data_file,
                          spectra_file=spectra_file,
                          layer=layers,
                          energy_min=energy_min,
                          energy_max=energy_max,
                          energy_step=energy_step,
                          repeat=repeat,
                          folder=folder,
                          baseline=baseline)

calibration.experiment.norm_to(norm_to_file)
calibration.experiment.slice(start=image_start, end=image_end)

calibrate_result = calibration.calibrate(source_to_detector_m=source_to_detector_m,
                                         offset_us=offset_us,
                                         vary='all',
                                         each_step=each_step)
calibration.index_peak(thres=0.1, min_dist=10)
# calibration.analyze_peak()
calibration.plot(before=False, index_level='iso', peak_id='indexed', peak_mark=False)
plt.show()
#
#
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
