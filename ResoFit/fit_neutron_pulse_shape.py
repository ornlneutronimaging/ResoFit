from ResoFit._pulse_shape import NeutronPulse
from lmfit import Parameters
from lmfit import minimize
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
from ResoFit.model import ikeda_carpenter
from ResoFit.model import cole_windsor
from ResoFit.model import cole_windsor_jparc
import pprint

# path1 = '/Users/Shawn/Dropbox (ORNL)/Postdoc_research/pulse_shape/source_section_1.dat'
# path2 = '/Users/Shawn/Dropbox (ORNL)/Postdoc_research/pulse_shape/source_section_2.dat'
path1 = '/Users/y9z/Dropbox (ORNL)/Postdoc_Research/pulse_shape/source_section_1.dat'
path2 = '/Users/y9z/Dropbox (ORNL)/Postdoc_Research/pulse_shape/source_section_2.dat'

neutron_pulse = NeutronPulse(path1)
neutron_pulse.load_shape_each(path2)
# neutron_pulse.export_total()
# neutron_pulse.export_each()
neutron_pulse.fit_shape(e_min=1, e_max=1000, model_index=1, drop=False, check_each=True, save_fig=True, save_df=False)
# pprint.pprint(neutron_pulse.shape_dict)

