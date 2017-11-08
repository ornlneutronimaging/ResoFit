from ResoFit._utilities import Layer
from ResoFit.calibration import Calibration
from ResoFit.fitresonance import FitResonance
import lmfit
import matplotlib.pyplot as plt
import pprint

# Global parameters
energy_min = 7
energy_max = 250
energy_step = 0.01
# Input sample name or names as str, case sensitive
# layer = 'UGd'
# thickness = 0.018  # mm
# density = get_foil_density_gcm3(length_mm=25, width_mm=25, thickness_mm=0.025, mass_g=0.14)
# density = None
# density = 8.86
layer_1 = 'U'
thickness_1 = 0.085
density_1 = None
layer_2 = 'Gd'
thickness_2 = 0.085
density_2 = None
# layer_3 = 'Cd'
# thickness_3 = 0.015
# density_3 = None
layer = Layer()
layer.add_layer(layer=layer_1, thickness_mm=thickness_1, density_gcm3=density_1)
layer.add_layer(layer=layer_2, thickness_mm=thickness_2, density_gcm3=density_2)
# layer.add_layer(layer=layer_3, thickness_mm=thickness_3, density_gcm3=density_3)

folder = 'data/IPTS_19558/reso_data_19558'
data_file = 'spheres.csv'
spectra_file = 'Image002_Spectra.txt'
image_start = None  # Can be omitted or =None
image_end = None  # Can be omitted or =None
norm_to_file = None  # 'sphere_background_1.csv'
baseline = True
each_step = False
before = False
table = True
grid = True
peak = 'indexed'
# items_to_plot = ['238-U', '235-U', 'Gd']
# items_to_plot = ['U-238', 'Gd-156', 'U']
# items_to_plot = [layer_1, layer_2]
# items_to_plot = [layer_2]
items_to_plot = None

repeat = 1
source_to_detector_m = 16.  # 16#16.445359069030175#16.447496101100739
offset_us = 0  # 0#2.7120797253959119#2.7355447625559037

# Calibrate source_to_detector and/or delay
calibration = Calibration(data_file=data_file,
                          spectra_file=spectra_file,
                          raw_layer=layer,
                          energy_min=energy_min,
                          energy_max=energy_max,
                          energy_step=energy_step,
                          repeat=repeat,
                          folder=folder,
                          baseline=baseline)

calibration.norm_to(norm_to_file)
calibration.slice(slice_start=300, slice_end=image_end)

calibrate_result = calibration.calibrate(source_to_detector_m=source_to_detector_m,
                                         offset_us=offset_us,
                                         vary='all',
                                         each_step=each_step)
# calibration.find_peak()
calibration.index_peak(thres=0.13, min_dist=21)

calibration.experiment.o_peak.analyze()

pprint.pprint(calibration.experiment.o_peak.peak_map_indexed)
pprint.pprint(calibration.experiment.o_peak.peak_df_scaled)


calibration.plot(before=before, table=table, peak=peak, grid=grid, items_to_plot=items_to_plot, interp=False)

# # Fit sample density or thickness
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
# # Fit isotope ratios
# # fit.fit_iso(layer=layer_1)
# fit.molar_conc()
# fit.index_peak(thres=0.15, min_dist=25)
# fit.plot(before=before, table=table, grid=grid, peak=peak,
#          items_to_plot=items_to_plot, interp=False)
#
# # fit.export()
