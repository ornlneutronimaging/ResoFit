import matplotlib.pyplot as plt

from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation

overwrite_csv = False
source_to_detector_m = 16.45

simulation = Simulation(energy_min=7, energy_max=150, energy_step=0.1, database='ENDF_VII')
simulation.add_layer(layer='Gd', layer_thickness_mm=0.075)
simulation._convolve_neutron_beam_shape(source_to_detector_m=source_to_detector_m,
                                        model_index=1,
                                        conv_proton=True,
                                        proton_params={})
# model_index:
# 1: 'ikeda_carpenter',
# 2: 'cole_windsor',
# 3: 'pseudo_voigt',
# 4: 'ikeda_carpenter_jparc',
# 5: 'cole_windsor_jparc'

folder = 'data/IPTS_19558/reso_data_19558'
data_file1 = 'Gd_thin.csv'
spectra_file = 'Image002_Spectra.txt'
experiment1 = Experiment(data_file=data_file1,
                         spectra_file=spectra_file,
                         folder=folder,
                         baseline=True)
# experiment1.slice(slice_start=300, reset_index=False)
# peak_df = experiment1.find_peak()
simulation.plot_simu(x_type='time', source_to_detector_m=source_to_detector_m, offset_us=2.67)
plt.plot(simulation.x_tof_us - 2.9, simulation.y_att, label='Simulated Data')
experiment1.plot_raw(x_type='time', time_unit='us')
# plt.title(title)


# path1 = '/Users/Shawn/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_1.dat'
# path2 = '/Users/Shawn/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_2.dat'
# # path1 = '/Users/y9z/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_1.dat'
# # path2 = '/Users/y9z/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_2.dat'
#
# source_to_detector_m = 16.45  # m
# neutron_pulse = NeutronPulse(path1, model_index=1)
# neutron_pulse.load_shape_each(path2)
#
# neutron_pulse.fit_shape(e_min=1, e_max=500, drop=False, norm=True, check_each=False, save_fig=False,
#                         overwrite_csv=False)
# neutron_pulse.fit_params(check_each=False, loglog_fit=True, overwrite_csv=False)
#
# neutron_pulse.plot_tof_shape_interp(neutron_pulse._energy_list_dropped,
#                                     source_to_detector_m=source_to_detector_m,
#                                     conv_proton=True,
#                                     for_sum=True)
