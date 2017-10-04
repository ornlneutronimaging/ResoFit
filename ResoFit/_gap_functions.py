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


def y_gap_for_fitting(params, exp_x_interp, exp_y_interp, layer, energy_min, energy_max, energy_step, each_step=False):
    parvals = params.valuesdict()
    layer_density_gcm3 = parvals['density_gcm3']
    layer_thickness_mm = parvals['thickness_mm']
    simulation = Simulation(energy_min=energy_min,
                            energy_max=energy_max,
                            energy_step=energy_step)
    simulation.add_layer(layer=layer, layer_thickness_mm=layer_thickness_mm, layer_density_gcm3=layer_density_gcm3)
    simu_x, simu_y = simulation.xy_simu(angstrom=False, transmission=False)
    gap = (exp_y_interp - simu_y)  # ** 2
    if each_step is True:
        print("density_gcm3: {}    thickness_mm: {}    chi^2: {}".format(layer_density_gcm3,
                                                                         layer_thickness_mm,
                                                                         sum((exp_y_interp - simu_y) ** 2)))
    return gap


def y_gap_for_fitting_multi_layers(params, exp_x_interp, exp_y_interp, energy_min, energy_max, energy_step, layers=[]):
    parvals = params.valuesdict()
    layer_density = parvals['density']
    layer_thickness = parvals['thickness']
    simulation = Simulation(energy_min=energy_min,
                            energy_max=energy_max,
                            energy_step=energy_step)
    for _i in range(len(layers)):
        simulation.add_layer(layer=layers[_i], layer_thickness_mm=layer_thickness, layer_density_gcm3=layer_density)
    simu_x, simu_y = simulation.xy_simu(angstrom=False, transmission=False)
    gap = (exp_y_interp - simu_y) ** 2
    return gap
