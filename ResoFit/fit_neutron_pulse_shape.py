from ResoFit._pulse_shape import NeutronPulse
from ResoFit._pulse_shape import ProtonPulse
from ResoFit.experiment import Experiment
import numpy as np
from ResoFit.simulation import Simulation
overwrite_csv = False
import scipy.signal as ss
import matplotlib.pyplot as plt

# path1 = '/Users/Shawn/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_1.dat'
# path2 = '/Users/Shawn/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_2.dat'
path1 = '/Users/y9z/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_1.dat'
path2 = '/Users/y9z/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_2.dat'

neutron_pulse = NeutronPulse(path1, model_index=1)

neutron_pulse.load_shape_each(path2)

# neutron_pulse.plot_total()

neutron_pulse.fit_shape(e_min=1, e_max=500, drop=False, norm=True, check_each=False, save_fig=False,
                        overwrite_csv=overwrite_csv)

neutron_pulse.fit_params(check_each=False, loglog_fit=True, overwrite_csv=overwrite_csv)

# e = np.linspace(5, 100, 20)
# neutron_pulse.make_shape([0.25, 1, 3, 4, 5, 6, 50, 100, 200])
# neutron_pulse.make_shape(e)

# neutron_pulse.plot_shape_total()


e_list = np.linspace(1, 100, 100)
# t_new = np.linspace(0.1, 30, 300)
t_new = None
# neutron_pulse.plot_shape_each_compare(e_min=15, e_max=200, norm=False, t_interp=None)
# neutron_pulse.plot_shape_interp(e_ev=e_list, t_interp=t_new, logy=False, norm=False)
neutron_pulse.plot_tof_shape_interp(e_ev=e_list, t_interp=t_new, for_sum=True, logy=False, norm=False)
# plt.show()
# neutron_pulse._make_shape(e_ev=e_list, t_interp=t_new, for_sum=True, norm=False, convolve_proton=True)
# neutron_pulse.shape_tof_df_interp.set_index('tof_us').sum(axis=1).plot()
# neutron_pulse._make_shape(e_ev=e_list, t_interp=t_new, for_sum=True, norm=False)


# simulation = Simulation(energy_min=7, energy_max=150, energy_step=0.1, database='ENDF_VII')
# simulation.add_layer(layer='Gd', layer_thickness_mm=0.15)
# simulation._convolve_beam_shape()
#
# folder = 'data/IPTS_19558/reso_data_19558'
# data_file1 = 'Gd_thick.csv'
# spectra_file = 'Image002_Spectra.txt'
# experiment1 = Experiment(data_file=data_file1,
#                          spectra_file=spectra_file,
#                          folder=folder,
#                          baseline=True)
# # experiment1.slice(slice_start=300, reset_index=False)
# # peak_df = experiment1.find_peak()
# simulation.plot_simu(x_type='time', source_to_detector_m=16.45, offset_us=2.67)
# plt.plot(simulation.x_tof_us-2.9, simulation.y_att, label='Convolution')
# experiment1.plot_raw(x_type='time', time_unit='us')
