from ResoFit.simulation import Simulation


def y_gap_for_calibration(params, simu_x, simu_y, energy_min, energy_max, energy_step, experiment,
                          baseline=False, each_step=False):
    # Unpack Parameters:
    parvals = params.valuesdict()
    source_to_detector_m = parvals['source_to_detector_m']
    offset_us = parvals['offset_us']
    exp_x, exp_y = experiment.xy_scaled(energy_min=energy_min,
                                        energy_max=energy_max,
                                        energy_step=energy_step,
                                        angstrom=False,
                                        transmission=False,
                                        offset_us=offset_us,
                                        source_to_detector_m=source_to_detector_m,
                                        baseline=baseline)

    gap = (exp_y - simu_y)  # ** 2
    if each_step is True:
        print("source_to_detector_m: {}    offset_us: {}    chi^2: {}".format(source_to_detector_m,
                                                                              offset_us,
                                                                              sum((exp_y - simu_y) ** 2)))
    return gap


def y_gap_for_fitting(params, exp_x_interp, exp_y_interp, layer_list,
                      energy_min, energy_max, energy_step,
                      each_step=False):
    parvals = params.valuesdict()
    simulation = Simulation(energy_min=energy_min,
                            energy_max=energy_max,
                            energy_step=energy_step)
    for each_layer in layer_list:
        simulation.add_layer(layer=each_layer,
                             layer_thickness_mm=parvals['thickness_mm_' + each_layer],
                             layer_density_gcm3=parvals['density_gcm3_' + each_layer])
    simu_x, simu_y = simulation.xy_simu(angstrom=False, transmission=False)
    gap = (exp_y_interp - simu_y)  # ** 2

    if each_step is True:
        for each_layer in layer_list:
            print("density_gcm3_{}: {}    thickness_mm_{}: {}    chi^2: {}".format(
                each_layer,
                parvals['density_gcm3_' + each_layer],
                each_layer,
                parvals['thickness_mm_' + each_layer],
                sum((exp_y_interp - simu_y) ** 2)))
    return gap


def y_gap_for_iso_fitting(params, exp_x_interp, exp_y_interp, layer, formatted_isotope_list,
                          fitted_simulation=Simulation(),
                          each_step=False):
    parvals = params.valuesdict()
    isotope_ratio_list = []
    for _isotope_index in range(len(formatted_isotope_list)):
        isotope_ratio_list.append(parvals[formatted_isotope_list[_isotope_index]])

    fitted_simulation.set_isotopic_ratio(layer=layer, element=layer, new_isotopic_ratio_list=isotope_ratio_list)
    simu_x, simu_y = fitted_simulation.xy_simu(angstrom=False, transmission=False)
    gap = (exp_y_interp - simu_y)  # ** 2
    return gap
