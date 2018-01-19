from ResoFit._pulse_shape import NeutronPulse
from lmfit import Parameters
from lmfit import minimize
from ResoFit._gap_functions import gap_neutron_pulse_ikeda_carpenter
from ResoFit._gap_functions import gap_neutron_pulse_cole_windsor
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
from ResoFit.model import ikeda_carpenter

# path1 = '/Users/Shawn/Dropbox (ORNL)/Postdoc_research/pulse_shape/source_section_1.dat'
# path2 = '/Users/Shawn/Dropbox (ORNL)/Postdoc_research/pulse_shape/source_section_2.dat'
path1 = '/Users/y9z/Dropbox (ORNL)/Postdoc_Research/pulse_shape/source_section_1.dat'
path2 = '/Users/y9z/Dropbox (ORNL)/Postdoc_Research/pulse_shape/source_section_2.dat'

neutron_pulse = NeutronPulse(path1)
neutron_pulse.load_shape_each(path2)
# neutron_pulse.export_total()
# neutron_pulse.export_each()
f = np.array(neutron_pulse.shape_dict[1]['f'])
t = np.array(neutron_pulse.shape_dict[1]['t_us'])

my_model = Model(ikeda_carpenter)
print("Parameters: ", my_model.param_names)
print("Independent variable: ", my_model.independent_vars)


my_model.set_param_hint('alpha', value=0.699, min=0, max=20)
my_model.set_param_hint('beta', value=0.0215, min=0, max=1)
my_model.set_param_hint('fraction', value=0.383, min=0, max=1)
my_model.set_param_hint('t0', value=0.0889, min=0, max=5)
my_model.set_param_hint('magnitude', value=1.46e12, min=0)
params = my_model.make_params(verbose=True)

my_model.print_param_hints()
print("+----------------- Fit neutron pulse shape -----------------+\nParams before:")
params.pretty_print()

# print(f)
# print(t)
each_step = True
# len(f)

# params_to_fitshape = Parameters()
# params_to_fitshape.add('alpha', value=0.06)
# params_to_fitshape.add('beta', value=0.05)
# params_to_fitshape.add('fraction', value=0.5, min=0, max=1)
# params_to_fitshape.add('t0', value=0.01, min=0)

# Print before
# print("+----------------- Fit neutron pulse shape -----------------+\nParams before:")
# params_to_fitshape.pretty_print()

# result = my_model.fit(f, params_to_fitshape, t=t, fit_kws={'nan_policy': 'omit'})

result = my_model.fit(f, params, t=t)


# shape_result = minimize(gap_neutron_pulse_ikeda_carpenter,
#                         params,
#                         method='leastsq',
#                         args=(t, f, each_step)
#                         )
