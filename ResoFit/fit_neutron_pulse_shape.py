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
neutron_pulse.fit_shape(e_min=1, e_max=1000, model_index=1, check_each=False, save_fig=False, save_df=True)
# pprint.pprint(neutron_pulse.shape_dict)





# f = np.array(neutron_pulse.shape_dict[1]['f_norm'])
# t = np.array(neutron_pulse.shape_dict[1]['t_us'])
#
# # my_model = Model(cole_windsor_jparc)
# my_model = Model(cole_windsor)
# print("Parameters: ", my_model.param_names)
# print("Independent variable: ", my_model.independent_vars)
#
# my_model.set_param_hint('sig1', value=0.06917, min=0, max=20)
# my_model.set_param_hint('sig2', value=0.2041, min=0, max=20)
# my_model.set_param_hint('gamma', value=6.291, min=0, max=20)
# # my_model.set_param_hint('gam1', value=6.291, min=0, max=20)
# # my_model.set_param_hint('gam2', value=1.285, min=0, max=20)
# my_model.set_param_hint('fraction', value=0.1308, min=0, max=1)
# my_model.set_param_hint('t0', value=0.3176, min=0, max=20)
# my_model.set_param_hint('norm_factor', value=0.9951, min=0)
#
# params = my_model.make_params(verbose=True)
#
# result = my_model.fit(f, params, t=t)
# result.params.pretty_print()
#
# print(result.fit_report())
#
# result.plot()
# plt.xlim(xmin=0, xmax=10)
# plt.show()

# my_model.set_param_hint('alpha', value=0.699, min=0, max=20)
# my_model.set_param_hint('beta', value=0.0215, min=0, max=20)
# my_model.set_param_hint('fraction', value=0.383, min=0, max=1)
# my_model.set_param_hint('t0', value=0.0889, min=0, max=20)
# # my_model.set_param_hint('norm_factor', value=1.46e12, min=0)
# # my_model.set_param_hint('norm_factor', value=1, min=0, vary=False)
# my_model.set_param_hint('norm_factor', value=1, min=0)
# params = my_model.make_params(verbose=True)
#
# # my_model.print_param_hints()
# print("+----------------- Fit neutron pulse shape -----------------+\nParams before:")
# params.pretty_print()

# print(f)
# print(t)
# each_step = True
# len(f

# params_to_fitshape = Parameters()
# params_to_fitshape.add('alpha', value=0.06)
# params_to_fitshape.add('beta', value=0.05)
# params_to_fitshape.add('fraction', value=0.5, min=0, max=1)
# params_to_fitshape.add('t0', value=0.01, min=0)

# Print before
# print("+----------------- Fit neutron pulse shape -----------------+\nParams before:")
# params_to_fitshape.pretty_print()

# result = my_model.fit(f, params_to_fitshape, t=t, fit_kws={'nan_policy': 'omit'})

# result = my_model.fit(f, params, t=t)
# print("+----------------- Fit results -----------------+\nParams before:")
# result.params.pretty_print()


# shape_result = minimize(gap_neutron_pulse_ikeda_carpenter,
#                         params,
#                         method='leastsq',
#                         args=(t, f, each_step)
#                         )
