# from scipy.optimize import leastsq

import pprint
from ImagingReso.resonance import Resonance

# Global parameters
_energy_min = 1e-5
_energy_max = 300
_energy_step = 0.01
# Input sample name or names as str, case sensitive
_layer_1 = 'Co'
_thickness_1 = 0.025 # mm
# _density_1 = 8 # g/cm3 deviated due to porosity

_layer_2 = 'Ag'
_thickness_2 = 0.025 # mm

o_reso = Resonance(energy_min=_energy_min, energy_max=_energy_max, energy_step=_energy_step)
o_reso.add_layer(formula=_layer_1, thickness=_thickness_1)
o_reso.add_layer(formula=_layer_2, thickness=_thickness_2)

# Ideal
pprint.pprint(o_reso.stack_sigma['Ag']['Ag']['energy_eV'])
pprint.pprint(o_reso.stack_sigma['Ag']['Ag']['sigma_b'])
pprint.pprint(o_reso.stack['Ag']['Ag']['isotopes']['list'])

energy_x = o_reso.stack_sigma['Ag']['Ag']['energy_eV']
sigma_1 = o_reso.stack_sigma['Ag']['Ag']['sigma_b']
sigma_2 = o_reso.stack_sigma['Co']['Co']['sigma_b']


# Get sigma dictionary


# '''Get atoms_per_cm^3 for each elements'''
# atoms_per_cm3_dict = {}
# for ele in elements:
#     atoms_per_cm3_dict[ele] = avogadro_number * density_gcm3_dict[ele]/molar_mass_dict[ele]
# print('Number of atoms per unit volume (#/cm^3) : ', atoms_per_cm3_dict)
#
# params = _fit_functions.def_params_from_dict(thick_cm_dict, 'thick_cm_')
# print(params)
# params = _fit_functions.add_params_from_doct(params, density_gcm3_dict, 'density_gcm3_')
# print(params)
# params = _fit_functions.add_params_from_doct(params, atoms_per_cm3_dict, 'atoms_per_cm3_')
# print(params)

# thick_cm_dict = _functions.dict_value_by_key(elements, thick_cm_list)
# print(elements)
# print(isotope_dict)
# print(sigma_iso_ele_eleisodict)
# print(sigma_iso_ele_sum_eledict)

# thick_cm_list = np.array(_input_thick_mm_list) / 10
# mixed_l_n_avo = _plot_functions.l_x_n_multi_ele_stack(elements,
#                                                       thick_cm_dict,
#                                                       density_gcm3_dict,
#                                                       molar_mass_dict)
#
# yi_values = list(dict.values(sigma_iso_ele_sum_eledict))
# yi_values_sum = sum(yi_values)
# # sum of (sigma * ele_ratio * iso_ratio)
# # print(yi_values)
# trans_sum = _functions.sig_l_2trans_quick(mixed_l_n_avo, yi_values_sum)
# y_trans_tot = trans_sum
#
# y_attenu_tot = 1 - y_trans_tot
#
# ideal_y_index = pku.indexes(y_attenu_tot, thres=0.15, min_dist=10)#, thres=0.1, min_dist=50)
# ideal_x_index = pku.interpolate(x_energy, y_attenu_tot, ind=ideal_y_index)
# print('x_ideal_peak: ', ideal_x_index)
# #
# plt.plot(x_energy, y_attenu_tot, 'b-', label=_input_ele_str+'_ideal')
# plt.plot(x_energy[ideal_y_index], y_attenu_tot[ideal_y_index], 'bo', label='peak_ideal')
#

# # Experiment
# source_to_detector_cm = 1610.1796635603498  # cm
# delay_ms = -12.112641168682274#-12.145 #4.5 - 16.61295379  # ms
# delay_us = delay_ms * 1000
# range_min = 500
# range_max = 2771
# _slice = range_min
# energy_min = 0
# time_lamda_ev_axis = 'eV'
# _name = 'foil6'
# data_path = 'data/' + _name + '.csv'
# spectra_path = 'data/spectra.txt'
# x_data_array = _functions.get_spectra_range(spectra_path, delay_us,
#                                             source_to_detector_cm, range_min, range_max)
# # print('x_exp: ', x_data_array)
# y_data_array = 1 - _functions.get_normalized_data_range(data_path, range_min, range_max)/4.25
# # print('y_exp: ', y_data_array)
# exp_y_index = pku.indexes(y_data_array, thres=0.12/max(y_data_array), min_dist=7)
# exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
# print('x_exp_peak: ', exp_x_index)
# print('Equal size: ', len(ideal_x_index) == len(exp_x_index))
# # baseline = pku.baseline(y_data_array)
# # print(baseline)
# # print('y_exp_peak: ', exp_y_index)
# #
# # df = pd.DataFrame()
# # df['Exp_x'] = x_data_array
# # df['Exp_y'] = y_data_array
# # df2 = pd.DataFrame()
# # df2['Ideal_x'] = x_energy
# # df2['Ideal_y'] = y_attenu_tot
# #
# # params = Parameters()
# # params.add('source_to_detector_cm', value=source_to_detector_cm)
# # params.add('delay_us', value=delay_us)
#
# # x_gap = _fit_funtions.peak_x_gap(params, ideal_x_index, y_data_array)
# # print('x_gap:', x_gap)
#
# # out = minimize(_fit_funtions.peak_x_gap, params, method='leastsq', args=(ideal_x_index, y_data_array))
# # out = scipy.optimize.minimize(_fit_funtions.peak_x_gap_scipy, delay_us, method='leastsq', args=(ideal_x_index, y_data_array))
# # print(out.residual)
#
# plt.plot(x_data_array, y_data_array, 'r-', label=_name)
# plt.plot(x_data_array[exp_y_index], y_data_array[exp_y_index], 'go', label='peak_exp')
# plt.title('Peak estimation')
#
# plt.ylim(-0.01, 1.01)
# plt.xlim(0, energy_max)
# plt.legend(loc='best')
# plt.show()
