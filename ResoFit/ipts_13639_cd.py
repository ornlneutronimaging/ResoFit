from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
from ResoFit.experiment import Experiment
import matplotlib.pyplot as plt
import numpy as np
import pprint
from ResoFit._utilities import get_foil_density_gcm3

# Global parameters
energy_min = 7
energy_max = 40
energy_step = 0.01
# Input sample name or names as str, case sensitive
layer_1 = 'Cd'
thickness_1 = 0.5  # mm
# density = get_foil_density_gcm3(length_mm=25, width_mm=25, thickness_mm=0.025, mass_g=0.14)
density = np.NaN
# density = 8.86

folder = 'data'
data_file = 'Cd.csv'
spectra_file = 'spectra.csv'
image_start = 500  # Can be omitted or =None
image_end = 2200  # Can be omitted or =None
norm_to_file = 'In.csv'
baseline = True
each_step = False

repeat = 1
source_to_detector_m = 16.123278721983177  # 16#16.445359069030175#16.447496101100739
offset_us = -12112.494119089204  # 0#2.7120797253959119#2.7355447625559037

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
                          folder=folder,
                          baseline=baseline)

calibration.norm_to(norm_to_file)
calibration.slice(slice_start=image_start, slice_end=image_end)

calibrate_result = calibration.calibrate(source_to_detector_m=source_to_detector_m,
                                         offset_us=offset_us,
                                         vary='all',
                                         each_step=each_step)
calibration.plot_before()
calibration.plot_after()

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
                   norm_to_file=norm_to_file,
                   slice_start=image_start,
                   slice_end=image_end,
                   baseline=baseline)
fit.fit(thickness_mm=thickness_1, density_gcm3=density, vary='thickness', each_step=each_step)
fit.molar_conc(layer_1)
fit.plot_before()
fit.plot_after()
