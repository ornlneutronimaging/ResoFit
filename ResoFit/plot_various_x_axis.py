from ImagingReso.resonance import Resonance

# Global parameters
_energy_min = 10
_energy_max = 170
_energy_step = 0.01
# Input sample name or names as str, case sensitive
_layer_1 = 'Gd'
_name = _layer_1
_thickness_1 = 0.075 # mm
# _density_1 = 8 # g/cm3 deviated due to porosity

o_reso = Resonance(energy_min=_energy_min, energy_max=_energy_max, energy_step=_energy_step)
o_reso.add_layer(formula=_layer_1, thickness=_thickness_1)

# o_reso.plot(all_elements=True, transmission=False, x_axis='time')
# o_reso.plot(all_elements=True, transmission=False, x_axis='lambda')
o_reso.plot(all_elements=True, transmission=False, x_axis='number')
