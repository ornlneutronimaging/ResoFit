from ImagingReso.resonance import Resonance
import numpy as np
# Global parameters
_energy_min = 13
_energy_max = 170
_energy_step = 0.01
# Input sample name or names as str, case sensitive
layer_1 = 'U'
thickness_1 = 0.05
density_1 = None
layer_2 = 'Gd'
thickness_2 = 0.05
density_2 = None

o_reso = Resonance(energy_min=_energy_min, energy_max=_energy_max, energy_step=_energy_step)
o_reso.add_layer(formula=layer_1, thickness=thickness_1, density=density_1)
o_reso.add_layer(formula=layer_2, thickness=thickness_2, density=density_2)

# o_reso.plot(all_elements=True, transmission=False, x_axis='time')
# o_reso.plot(all_elements=True, transmission=False, x_axis='lambda')
o_reso.plot(mixed=True, all_isotopes=False, y_axis='attenuation', x_axis='energy', offset_us=0, time_resolution_us=0.16,
            source_to_detector_m=16.125)
