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
energy_min = 4.1
energy_max = 800
energy_step = 0.01
# Input sample name or names as str, case sensitive
layers = Layer()
layers.add_layer(layer='Ta', thickness_mm=0.127)

folder = 'data/IPTS_20440'
spectra_file = 'spectra.txt'
data_file = 'Ta_80C_12pC.csv'
image_start = None  # Can be omitted or =None
image_end = None  # Can be omitted or =None
norm_to_file = 'OB_80C_12pC.csv'

baseline = False
each_step = False

norm_factor = 1.
source_to_detector_m = 16.45  # 16#16.445359069030175#16.447496101100739
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
                                         vary='source_to_detector',
                                         # vary='all',
                                         each_step=each_step)
calibration.index_peak(thres=0.05, min_dist=2, map_min_dist=5, map_thres=0.05)
# calibration.analyze_peak()
calibration.experiment.plot()
calibration.plot(y_type='attenuation',
                 # y_type='transmission',
                 x_type='energy',
                 # t_unit='ms',
                 before=True,
                 # interp=True,
                 mixed=True,
                 # peak_exp='all',
                 table=False,
                 peak_exp='indexed',
                 peak_height=True,
                 index_level='ele',
                 peak_id='all',
                 logx=False,
                 )
plt.xlim(left=0, right=300)
plt.show()

# # Fit the peak height
# fit = FitResonance(spectra_file=spectra_file,
#                    data_file=data_file,
#                    folder=folder,
#                    repeat=norm_factor,
#                    energy_min=energy_min,
#                    energy_max=energy_max,
#                    energy_step=energy_step,
#                    calibrated_offset_us=calibration.calibrated_offset_us,
#                    calibrated_source_to_detector_m=calibration.calibrated_source_to_detector_m,
#                    norm_to_file=norm_to_file,
#                    slice_start=image_start,
#                    slice_end=image_end,
#                    baseline=baseline)
# fit_result = fit.fit('Ta', vary='density', each_step=each_step)
# fit.molar_conc()
# fit.index_peak(thres=0.15, min_dist=25)
# # fit.fit_iso(layer=layer_2)
# fit.plot(peak_id='all', interp=False)
# # fit.export('Exp_Gd_150_um.csv')
