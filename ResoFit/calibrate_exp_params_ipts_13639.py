from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
from ResoFit.experiment import Experiment
import matplotlib.pyplot as plt
import numpy as np
import pprint

# Global parameters
energy_min = 7
energy_max = 150
energy_step = 0.01
# Input sample name or names as str, case sensitive
layer_1 = 'Gd'
thickness_1 = 0.15  # mm
mass = 0.36  # gram
length = 25
width = 25
height = 0.075
mm3_to_cm3 = 0.001
density = np.NaN  # mass / (length * width * height * mm3_to_cm3)

folder = 'data'
data_file = 'Ag.csv'
spectra_file = 'spectra.csv'

repeat = 5
source_to_detector_m = 16.12  # 16#16.445359069030175#16.447496101100739
offset_us = 12  # 0#2.7120797253959119#2.7355447625559037

# Calibrate the peak positions
experiment = Experiment(data_file=data_file,
                        spectra_file=spectra_file,
                        repeat=repeat,
                        folder=folder)
# exp_x, exp_y = experiment.xy_scaled(energy_min=energy_min,
#                                     energy_max=energy_max,
#                                     energy_step=energy_step,
#                                     offset_us=offset_us,
#                                     source_to_detector_m=source_to_detector_m)

exp_x_sliced, exp_y_sliced = experiment.xy_sliced(450, 2000, baseline=False)
plt.plot(exp_y_sliced)
plt.show()
