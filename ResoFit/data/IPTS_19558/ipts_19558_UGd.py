from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
from ResoFit.experiment import Experiment
import matplotlib.pyplot as plt
import numpy as np
import pprint
from ResoFit._utilities import get_foil_density_gcm3
from ResoFit._utilities import Layer

# Global parameters
energy_min = 7
energy_max = 300
energy_step = 0.01

layers = Layer()
layers.add_layer(layer='U', thickness_mm=0.018, density_gcm3=None)
layers.add_layer(layer='Gd', thickness_mm=0.015, density_gcm3=None)

folder = 'data/IPTS_19558/reso_data_19558'
data_file = 'spheres.csv'
spectra_file = 'Image002_Spectra.txt'
database = 'ENDF_VIII'
image_start = 300  # Can be omitted or =None
image_end = None  # Can be omitted or =None
norm_to_file = None  # 'sphere_background_1.csv'
# norm_to_file = 'sphere_background_1.csv'
baseline = True
each_step = False
baseline_deg = 3

norm_factor = 1
source_to_detector_m = 16.43  # 16#16.445359069030175#16.447496101100739
offset_us = 2.7  # 0#2.7120797253959119#2.7355447625559037

# x_type = 'lambda'
# x_type = 'energy'
x_type = 'number'
# x_type = 'time'
# y_type = 'transmission'
y_type = 'attenuation'

# Calibrate the peak positions
calibration = Calibration(data_file=data_file,
                          spectra_file=spectra_file,
                          layer=layers,
                          energy_min=energy_min,
                          energy_max=energy_max,
                          energy_step=energy_step,
                          folder=folder,
                          exp_source_to_detector_m=source_to_detector_m,
                          exp_offset_us=offset_us,
                          database=database,
                          baseline=baseline,
                          baseline_deg=baseline_deg,
                          x_type=x_type,
                          y_type=y_type
                          )

calibration.experiment.norm_to(file=norm_to_file, norm_factor=norm_factor)
calibration.experiment.slice(start=image_start, end=image_end)

calibrate_result = calibration.calibrate(source_to_detector_m=source_to_detector_m,
                                         offset_us=offset_us,
                                         vary='all',
                                         each_step=each_step)
calibration.index_peak(thres_exp=0.15, min_dist_exp=20, min_dist_map=15, thres_map=0.12, rel_tol=0.01)
calibration.analyze_peak(report=False, fit_model='Lorentzian')  # ['Gaussian', 'Lorentzian']


calibration.plot(y_type=y_type,
                 x_type=x_type,
                 # t_unit='ns',
                 # before=True,
                 # interp=True,
                 mixed=True,
                 table=True,
                 peak_exp='all',
                 peak_height=True,
                 index_level='ele',
                 # peak_id='all',
                 logx=False,
                 )
# plt.xlim(left=0, right=400)
plt.show()
# calibration.plot(before=before, items_to_plot=items_to_plot)
# calibration.plot(before=before, items_to_plot=['Gd', 'U*', '235-U'])

# # Fit the peak height
# fit = FitResonance(spectra_file=spectra_file,
#                    data_file=data_file,
#                    folder=folder,
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
# fit_result = fit.fit(layer, vary='density', each_step=each_step)
# fit.molar_conc()
# # fit.fit_iso(layer=layer_1)
# # fit.molar_conc()
# fit.plot(before=before, items_to_plot=['Gd', 'U*', '235-U'])
