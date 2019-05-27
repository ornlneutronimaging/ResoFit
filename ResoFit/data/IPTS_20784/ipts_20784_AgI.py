from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
from ResoFit.experiment import Experiment
import matplotlib.pyplot as plt
import numpy as np
import pprint
from ResoFit._utilities import get_foil_density_gcm3
from ResoFit._utilities import Layer
import lmfit

# Global parameters
energy_min = 3
energy_max = 200
energy_step = 0.01
# Input sample name or names as str, case sensitive
layers = Layer()
layers.add_layer(layer='Ag', thickness_mm=0.0635)
layers.add_layer(layer='I', thickness_mm=0.0635)

folder = 'data/IPTS_20784/reso_data_20784'
# data_file2 = 'spheres_background_1.csv'
spectra_file = 'Ta_lead_10mm__0__040_Spectra.txt'
# data_file = 'AgI.csv'
data_file = 'AgI_pellets_all.csv'
image_start = None  # Can be omitted or =None
image_end = None  # Can be omitted or =None
# norm_to_file = 'blank_region.csv'
norm_to_file = 'blank_pellets_all.csv'
baseline = True
each_step = False

norm_factor = 1
source_to_detector_m = 16.5  # 16#16.445359069030175#16.447496101100739
offset_us = 0  # 0#2.7120797253959119#2.7355447625559037

# Calibrate the peak positions
calibration = Calibration(data_file=data_file,
                          spectra_file=spectra_file,
                          layer=layers,
                          energy_min=energy_min,
                          energy_max=energy_max,
                          energy_step=energy_step,
                          folder=folder,
                          baseline=baseline)

calibration.experiment.norm_to(norm_to_file, norm_factor=norm_factor)
calibration.experiment.slice(start=image_start, end=image_end)
calibrate_result = calibration.calibrate(source_to_detector_m=source_to_detector_m,
                                         offset_us=offset_us,
                                         vary='all',
                                         # vary='source_to_detector',
                                         each_step=each_step)
calibration.index_peak(thres=0.05, min_dist=2, map_min_dist=5, map_thres=0.05)
# calibration.analyze_peak()
calibration.experiment.plot()
calibration.plot(y_type='attenuation',
                 # y_type='transmission',
                 x_type='energy',
                 # t_unit='ms',
                 # before=True,
                 # interp=True,
                 # mixed=True,
                 # peak_exp='all',
                 table=False,
                 # peak_exp='indexed',
                 peak_height=False,
                 index_level='ele',
                 peak_id='indexed',
                 logx=False,
                 )
plt.xlim(left=1, right=100)
plt.show()

calibration.export(y_type='attenuation',
                   # y_type='transmission',
                   x_type='energy',
                   # t_unit='ms',
                   # before=True,
                   # interp=True,
                   # mixed=True,
                   # peak_exp='all',
                   # peak_exp='indexed',
                   index_level='ele',
                   peak_id='indexed',
                   )
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
# fit_result = fit.fit(layer, vary=fit_vary, each_step=each_step)
# fit.molar_conc()
# fit.index_peak(thres=0.15, min_dist=25)
# # fit.fit_iso(layer=layer_2)
# fit.plot(peak_id='all', interp=False)
# # fit.export('Exp_Gd_150_um.csv')