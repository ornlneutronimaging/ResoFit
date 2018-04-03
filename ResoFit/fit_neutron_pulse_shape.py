from ResoFit._pulse_shape import NeutronPulse
from ResoFit._pulse_shape import ProtonPulse
from ResoFit.experiment import Experiment
import numpy as np
from ResoFit.simulation import Simulation
overwrite_csv = False
import scipy.signal as ss
import matplotlib.pyplot as plt


source_to_detector_m = 16.45

simulation = Simulation(energy_min=7, energy_max=150, energy_step=0.1, database='ENDF_VII')
simulation.add_layer(layer='Gd', layer_thickness_mm=0.15)
simulation._convolve_neutron_beam_shape(source_to_detector_m=source_to_detector_m)

folder = 'data/IPTS_19558/reso_data_19558'
data_file1 = 'Gd_thick.csv'
spectra_file = 'Image002_Spectra.txt'
experiment1 = Experiment(data_file=data_file1,
                         spectra_file=spectra_file,
                         folder=folder,
                         baseline=True)
# experiment1.slice(slice_start=300, reset_index=False)
# peak_df = experiment1.find_peak()
simulation.plot_simu(x_type='time', source_to_detector_m=source_to_detector_m, offset_us=2.67)
plt.plot(simulation.x_tof_us-2.9, simulation.y_att, label='Convolution')
experiment1.plot_raw(x_type='time', time_unit='us')
# plt.title(title)
plt.show()


# proton_path = '/Users/y9z/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/proton_pulse/waveform_20170901.txt'
#
# proton_pulse = ProtonPulse(path=proton_path)
