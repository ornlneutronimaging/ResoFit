from ResoFit.simulation import Simulation
import ImagingReso._utilities as reso_util
from ResoFit.model import ikeda_carpenter
from ResoFit.model import cole_windsor


def y_gap_for_calibration(params, simu_x, simu_y, energy_min, energy_max, energy_step, experiment,
                          x_type, y_type, baseline=False, each_step=False):
    # Unpack Parameters:
    parvals = params.valuesdict()
    source_to_detector_m = parvals['source_to_detector_m']
    offset_us = parvals['offset_us']
    exp_x, exp_y = experiment.xy_scaled(energy_min=energy_min,
                                        energy_max=energy_max,
                                        energy_step=energy_step,
                                        x_type=x_type,
                                        y_type=y_type,
                                        offset_us=offset_us,
                                        source_to_detector_m=source_to_detector_m,
                                        baseline=baseline)

    gap = abs(exp_y - simu_y)  # ** 2
    if each_step is True:
        print("Trying: source_to_detector_m: {}    offset_us: {}    chi^2: {}".format(source_to_detector_m,
                                                                                      offset_us,
                                                                                      sum((exp_y - simu_y) ** 2)))
    return gap


# def y_gap_for_adv_calibration(params, ideal_x, thres, min_dist, experiment, each_step=False):
#     # Unpack Parameters:
#     parvals = params.valuesdict()
#     source_to_detector_m = parvals['source_to_detector_m']
#     offset_us = parvals['offset_us']
#     exp_peak_df = experiment.find_peak(thres=thres, min_dist=min_dist)
#     exp_x = reso_util.s_to_ev(array=exp_peak_df['x_s'],
#                               source_to_detector_m=source_to_detector_m,
#                               offset_us=offset_us)
#     sorted(exp_x)
#     gap = (exp_x - ideal_x)  # ** 2
#     if each_step is True:
#         print("Trying: source_to_detector_m: {}    offset_us: {}    chi^2: {}".format(source_to_detector_m,
#                                                                                       offset_us,
#                                                                                       sum((exp_x - ideal_x) ** 2)))
#     return gap


def y_gap_for_fitting(params, exp_x_interp, exp_y_interp, layer_list,
                      energy_min, energy_max, energy_step, database,
                      each_step=False):
    parvals = params.valuesdict()
    simulation = Simulation(energy_min=energy_min,
                            energy_max=energy_max,
                            energy_step=energy_step,
                            database=database)
    for each_layer in layer_list:
        simulation.add_layer(layer=each_layer,
                             thickness_mm=parvals['thickness_mm_' + each_layer],
                             density_gcm3=parvals['density_gcm3_' + each_layer])
    simu_x = simulation.get_x(x_type='energy')
    simu_y = simulation.get_y(y_type='attenuation')
    gap = (exp_y_interp - simu_y)  # ** 2

    if each_step is True:
        for each_layer in layer_list:
            print("Trying: density_gcm3_{}: {}    thickness_mm_{}: {}    chi^2: {}".format(
                each_layer,
                parvals['density_gcm3_' + each_layer],
                each_layer,
                parvals['thickness_mm_' + each_layer],
                sum((exp_y_interp - simu_y) ** 2)))
    return gap


def y_gap_for_iso_fitting(params, exp_x_interp, exp_y_interp, layer, formatted_isotope_list,
                          fitted_simulation: Simulation,
                          each_step=False):
    parvals = params.valuesdict()
    isotope_ratio_list = []
    for _isotope_index in range(len(formatted_isotope_list)):
        isotope_ratio_list.append(parvals[formatted_isotope_list[_isotope_index]])

    fitted_simulation.set_isotopic_ratio(layer=layer, element=layer, new_isotopic_ratio_list=isotope_ratio_list)
    simu_x = fitted_simulation.get_x(x_type='energy')
    simu_y = fitted_simulation.get_y(y_type='attenuation')
    gap = (exp_y_interp - simu_y)  # ** 2

    if each_step is True:
        for each_iso in formatted_isotope_list:
            print("Trying: {}: {}    chi^2: {}".format(
                each_iso,
                parvals[each_iso],
                sum((exp_y_interp - simu_y) ** 2))
            )
    return gap

# def gap_neutron_pulse_ikeda_carpenter(params, t, f, each_step=False):
#     parvals = params.valuesdict()
#     alpha = parvals['alpha']
#     beta = parvals['beta']
#     fraction = parvals['fraction']
#     t0 = parvals['t0']
#     simulated_shape = ikeda_carpenter(t=t,
#                                       alpha=alpha,
#                                       beta=beta,
#                                       fraction=fraction,
#                                       t0=t0)
#     # print(simulated_shape)
#     gap = f - simulated_shape
#     # print(gap)
#
#     if each_step is True:
#         print("Trying: alpha: {}    beta: {}    fraction: {}    t0: {}     chi^2: {}".format(
#             parvals['alpha'], parvals['beta'], parvals['fraction'], parvals['t0'],
#             sum((f - simulated_shape) ** 2))
#         )
#     return gap
#
#
# def gap_neutron_pulse_cole_windsor(params, t, f, each_step=False):
#     parvals = params.valuesdict()
#     sig1 = parvals['sig1']
#     sig2 = parvals['sig2']
#     gam1 = parvals['gam1']
#     gam2 = parvals['gam2']
#     norm_factor = parvals['norm_factor']
#     fraction = parvals['fraction']
#     t0 = parvals['t0']
#     simulated_shape = cole_windsor(t=t,
#                                    sig1=sig1,
#                                    sig2=sig2,
#                                    gam1=gam1,
#                                    gam2=gam2,
#                                    norm_factor=norm_factor,
#                                    fraction=fraction,
#                                    t0=t0)
#     gap = f - simulated_shape
#
#     if each_step is True:
#         print(
#             "Trying: sig1: {}    sig2: {}    gam1: {}    gam2: {}    norm_factor: {}    fraction: {}    t0: {}     chi^2: {}".format(
#                 parvals['sig1'], parvals['sig2'],
#                 parvals['gam1'], parvals['gam2'],
#                 parvals['norm_factor'],
#                 parvals['fraction'],
#                 parvals['t0'],
#                 sum((f - simulated_shape) ** 2))
#         )
#     return gap
