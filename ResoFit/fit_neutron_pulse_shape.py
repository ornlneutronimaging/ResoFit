from ResoFit._pulse_shape import NeutronPulse


path1 = '/Users/Shawn/Dropbox (ORNL)/Postdoc_research/pulse_shape/source_section_1.dat'
path2 = '/Users/Shawn/Dropbox (ORNL)/Postdoc_research/pulse_shape/source_section_2.dat'

neutron_pulse = NeutronPulse(path1)
neutron_pulse.load_shape_each(path2)
