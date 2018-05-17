from ResoFit.simulation import Simulation
import pprint
import numpy as np

# Global parameters
energy_min = 7
energy_max = 150
energy_step = 0.01
# Input sample name or names as str, case sensitive
layer_1 = 'Gd'
thickness_1 = 0.15
# density_1 = None
layer_2 = 'Pb'
thickness_2 = 2.54
# density_2 = None

simulation = Simulation(energy_min=energy_min,
                        energy_max=energy_max,
                        energy_step=energy_step,
                        database='ENDF_VIII')

simulation.o_reso.add_Layer(formula=layer_1, thickness=thickness_1)
# simulation.o_reso.add_Layer(formula=layer_2, thickness=thickness_2)

# o_reso.plot(all_elements=True, transmission=False, x_axis='time')
# o_reso.plot(all_elements=True, transmission=False, x_axis='lambda')
pprint.pprint(simulation.o_reso.stack)
simulation.o_reso.plot(mixed=True,
                       all_elements=True,
                       all_isotopes=False,
                       y_axis='attenuation',
                       x_axis='energy', offset_us=0,
                       time_resolution_us=0.16,
                       source_to_detector_m=16.125)
