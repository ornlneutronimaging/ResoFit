from ResoFit._pulse_shape import NeutronPulse


# path1 = '/Users/Shawn/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_1.dat'
# path2 = '/Users/Shawn/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_2.dat'
path1 = '/Users/y9z/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_1.dat'
path2 = '/Users/y9z/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_2.dat'

neutron_pulse = NeutronPulse(path1, model_index=1)
neutron_pulse.load_shape_each(path2)

neutron_pulse.plot_total()

# neutron_pulse.fit_shape(e_min=1, e_max=500,
#                         drop=False,
#                         norm=True,
#                         check_each=False,
#                         save_fig=False,
#                         overwrite_csv=False)
#
# neutron_pulse.fit_params(check_each=True,loglog_fit=True)

