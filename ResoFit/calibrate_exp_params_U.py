from lmfit import Parameters
from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
import numpy as np
import pprint

# Global parameters
energy_min = 7
energy_max = 150
energy_step = 0.01
# Input sample name or names as str, case sensitive
_layer_1 = 'U'
_thickness_1 = 0.15  # mm
mass = 0.36  # gram
length = 25
width = 25
height = 0.075
mm3_to_cm3 = 0.001
density = mass / (length * width * height * mm3_to_cm3)

_density_1 = np.NaN
folder = 'data'
data_file = 'values_no_xy.txt'
spectra_file = 'Image002_Spectra.txt'

repeat = 1
source_to_detector_m = 16.  # 16#16.445359069030175#16.447496101100739
offset_us = 0  # 0#2.7120797253959119#2.7355447625559037

# Calibrate the peak positions
calibration = Calibration(data_file=data_file,
                          spectra_file=spectra_file,
                          layer_1=_layer_1,
                          thickness_1=_thickness_1,
                          density_1=np.NaN,
                          energy_min=energy_min,
                          energy_max=energy_max,
                          energy_step=energy_step,
                          repeat=repeat,
                          folder=folder)

calibrate_result = calibration.calibrate(source_to_detector_m=source_to_detector_m, offset_us=offset_us)
calibration.plot_before()
calibration.plot_after()
# calibration.plot_after_interp()

# Fit the peak height
fit = FitResonance(spectra_file=spectra_file, data_file=data_file, repeat=repeat, layer=_layer_1,
                   calibrated_offset_us=calibration.calibrated_offset_us,
                   calibrated_source_to_detector_m=calibration.calibrated_source_to_detector_m,
                   energy_min=energy_min, energy_max=energy_max, energy_step=energy_step)
fit.fit(thickness=_thickness_1, density=density, vary='density')
fit.molar_conc(_layer_1)
fit.plot_before()
fit.plot_after()
