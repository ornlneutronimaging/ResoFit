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

neutron_pulse.plot_shape_each_compare(e_min=15, e_max=100, norm=False)


