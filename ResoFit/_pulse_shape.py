import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import Parameters
from lmfit import minimize
from ResoFit._gap_functions import gap_neutron_pulse_ikeda_carpenter
from ResoFit._gap_functions import gap_neutron_pulse_cole_windsor


class NeutronPulse(object):

    def __init__(self, model_index=1):
        """

        :param model_index: [1: 'ikeda_carpenter', 2: 'cole_windsor', 3: 'pseudo_voigt']
        :type model_index: int
        """
        self.model_index = model_index
        self.params_to_fitshape = None
        self.shape_result = None

    def fit_shape(self, t, f, each_step=False):

        self.params_to_fitshape = Parameters()

        # ikeda_carpenter
        if self.model_index == 1:
            # Load params
            self.params_to_fitshape.add('alpha',
                                        # value=source_to_detector_m
                                        )
            self.params_to_fitshape.add('beta',
                                        # value=offset_us
                                        )
            self.params_to_fitshape.add('fraction',
                                        # value=0.5,
                                        min=0,
                                        max=1
                                        )
            self.params_to_fitshape.add('t0',
                                        # value=offset_us
                                        )
            # Use lmfit to obtain params by minimizing gap_function
            self.shape_result = minimize(gap_neutron_pulse_ikeda_carpenter,
                                         self.params_to_fitshape,
                                         method='leastsq',
                                         args=(t, f, each_step)
                                         )
        # cole_windsor
        elif self.model_index == 2:
            # Load params
            self.params_to_fitshape.add('sig1',
                                        # value=source_to_detector_m
                                        )
            self.params_to_fitshape.add('sig2',
                                        # value=offset_us
                                        )
            self.params_to_fitshape.add('gam1',
                                        # value=offset_us
                                        )
            self.params_to_fitshape.add('gam2',
                                        # value=offset_us
                                        )
            self.params_to_fitshape.add('norm_factor',
                                        # value=source_to_detector_m
                                        )
            self.params_to_fitshape.add('fraction',
                                        # value=0.5,
                                        min=0,
                                        max=1
                                        )
            self.params_to_fitshape.add('t0',
                                        # value=offset_us,
                                        vary=True
                                        )
            # Use lmfit to obtain params by minimizing gap_function
            self.shape_result = minimize(gap_neutron_pulse_cole_windsor,
                                         self.params_to_fitshape,
                                         method='leastsq',
                                         args=(t, f, each_step))

        # Print before
        print("+----------------- Fit neutron pulse shape -----------------+\nParams before:")
        self.params_to_fitshape.pretty_print()
        # Use lmfit to obtain params by minimizing gap_function

        # Print after
        print("\nParams after:")
        self.shape_result.__dict__['params'].pretty_print()
        # Print chi^2
        print("Calibration chi^2 : {}\n".format(self.shape_result.__dict__['chisqr']))


class ProtonPulse(object):
    pass
