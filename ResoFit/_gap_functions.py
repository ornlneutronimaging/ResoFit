from ResoFit.simulation import Simulation


def y_gap_for_calibration(params, simu_x, simu_y, energy_min, energy_max, energy_step, experiment, baseline=False):
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
    # print(sum((simu_x - exp_x) ** 2))
    # if sum((simu_x - exp_x) ** 2) > 1e-10:
    #     raise ValueError("The interpolated experiment x-axis is not identical to simulation x-axis!")
    gap = (exp_y - simu_y) #** 2
    return gap


def y_gap_for_fitting(params, exp_x_interp, exp_y_interp, layer, energy_min, energy_max, energy_step):
    parvals = params.valuesdict()
    layer_density = parvals['density']
    layer_thickness = parvals['thickness']
    simulation = Simulation(energy_min=energy_min,
                            energy_max=energy_max,
                            energy_step=energy_step)
    simulation.add_layer(layer=layer, layer_thickness=layer_thickness, layer_density=layer_density)
    simu_x, simu_y = simulation.xy_simu(angstrom=False, transmission=False)
    # if sum((simu_x - exp_x_interp) ** 2) > 1e-10:
    #     raise ValueError("The interpolated experiment x-axis is not identical to simulation x-axis!")
    gap = (exp_y_interp - simu_y) #** 2
    return gap


def y_gap_for_fitting_multi_layers(params, exp_x_interp, exp_y_interp, energy_min, energy_max, energy_step, layers=[]):
    parvals = params.valuesdict()
    layer_density = parvals['density']
    layer_thickness = parvals['thickness']
    simulation = Simulation(energy_min=energy_min,
                            energy_max=energy_max,
                            energy_step=energy_step)
    for _i in range(len(layers)):
        simulation.add_layer(layer=layers[_i], layer_thickness=layer_thickness, layer_density=layer_density)
    simu_x, simu_y = simulation.xy_simu(angstrom=False, transmission=False)
    gap = (exp_y_interp - simu_y) ** 2
    return gap
