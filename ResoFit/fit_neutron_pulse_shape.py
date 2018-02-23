from ResoFit._pulse_shape import NeutronPulse
import numpy as np

overwrite_csv = False
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


e_list = np.linspace(1, 300, 300)
t_new = np.linspace(0.1, 30, 300)
# neutron_pulse.plot_shape_each_compare(e_min=15, e_max=200, norm=False, t_interp=None)
# neutron_pulse.plot_tof_shape_interp(e_ev=e_list, t_interp=t_new, for_sum=True, logy=False, norm=False)
neutron_pulse._make_shape(e_ev=e_list, t_interp=t_new, for_sum=True, norm=False)
neutron_pulse.shape_tof_df_interp.set_index('tof_us').sum(axis=1).plot()
# neutron_pulse._make_shape(e_ev=e_list, t_interp=t_new, for_sum=True, norm=False)
