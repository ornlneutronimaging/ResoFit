from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
from ResoFit.experiment import Experiment
import matplotlib.pyplot as plt
import numpy as np
import pprint
from ResoFit._utilities import get_foil_density_gcm3
from ResoFit._utilities import Layer


# Global parameters
energy_min = 13
energy_max = 150
energy_step = 0.01
# Input sample name or names as str, case sensitive
# layer = 'UGd'
# thickness = 0.018  # mm
# density = get_foil_density_gcm3(length_mm=25, width_mm=25, thickness_mm=0.025, mass_g=0.14)
# density = None
# density = 8.86
layer_1 = 'U'
thickness_1 = 0.018
density_1 = None
# layer_2 = 'Gd'
# thickness_2 = 0.015
# density_2 = None
# layer_3 = 'Cd'
# thickness_3 = 0.015
# density_3 = None
layer = Layer()
layer.add_layer(layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
# layer.add_Layer(layer=layer_2, thickness_mm=thickness_2, density_gcm3=density_2)
# layer.add_Layer(layer=layer_3, thickness_mm=thickness_3, density_gcm3=density_3)

folder = 'data/IPTS_18521/reso_data_18521'
data_file = 'run_33_resonance.txt'
spectra_file = 'Image033_Spectra.txt'
image_start = None  # Can be omitted or =None
image_end = None  # Can be omitted or =None
norm_to_file = 'run_33_resonance_ob.txt'#None  # 'sphere_background_1.csv'
baseline = True
each_step = False

repeat = 1
source_to_detector_m = 15.4318  # 16#16.445359069030175#16.447496101100739
offset_us = 2.7  # 0#2.7120797253959119#2.7355447625559037

# Calibrate the peak positions
calibration = Calibration(data_file=data_file,
                          spectra_file=spectra_file,
                          layer=layer,
                          energy_min=energy_min,
                          energy_max=energy_max,
                          energy_step=energy_step,
                          norm_factor=repeat,
                          folder=folder,
                          baseline=baseline)

calibration.norm_to(norm_to_file)
calibration.slice(slice_start=image_start, slice_end=image_end)

calibrate_result = calibration.calibrate(source_to_detector_m=source_to_detector_m,
                                         offset_us=offset_us,
                                         vary='all',
                                         each_step=each_step)
calibration.plot()

# Fit the peak height
fit = FitResonance(spectra_file=spectra_file,
                   data_file=data_file,
                   folder=folder,
                   repeat=repeat,
                   energy_min=energy_min,
                   energy_max=energy_max,
                   energy_step=energy_step,
                   calibrated_offset_us=calibration.calibrated_offset_us,
                   calibrated_source_to_detector_m=calibration.calibrated_source_to_detector_m,
                   norm_to_file=norm_to_file,
                   slice_start=image_start,
                   slice_end=image_end,
                   baseline=baseline)
fit.fit(layer, vary='density', each_step=each_step)
fit.molar_conc()
fit.plot(error=True)
