import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import peakutils as pku
from lmfit import Parameters
from scipy.interpolate import interp1d
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation
import numpy as np
from lmfit import minimize
from ResoFit._gap_functions import y_gap_for_calibration


class Calibration(Simulation):
    def __init__(self, spectra_file, data_file, raw_layer,
                 energy_min=1e-5, energy_max=1000, energy_step=0.01,
                 repeat=1, folder='data', baseline=False):
        """
        Initialization with passed file location and sample info
        :param spectra_file:
        :type spectra_file:
        :param data_file:
        :type data_file:
        :param raw_layer:
        :type raw_layer:
        :param energy_min:
        :type energy_min:
        :param energy_max:
        :type energy_max:
        :param energy_step:
        :type energy_step:
        :param repeat:
        :type repeat:
        :param folder:
        :type folder:
        :param baseline: True -> to remove baseline/background by detrend
        :type baseline: boolean
        """
        super().__init__(energy_min=energy_min, energy_max=energy_max, energy_step=energy_step)
        for _each_layer in list(raw_layer.info.keys()):
            self.add_layer(layer=_each_layer,
                           layer_thickness_mm=raw_layer.info[_each_layer]['thickness']['value'],
                           layer_density_gcm3=raw_layer.info[_each_layer]['density']['value'])
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.energy_step = energy_step
        self.experiment = Experiment(spectra_file=spectra_file, data_file=data_file, repeat=repeat, folder=folder)
        self.repeat = repeat
        self.data_file = data_file
        self.spectra_file = spectra_file
        self.init_source_to_detector_m = None
        self.init_offset_us = None
        self.calibrated_offset_us = None
        self.calibrated_source_to_detector_m = None
        self.calibrate_result = None
        self.exp_x_raw_calibrated = None
        self.exp_y_raw_calibrated = None
        self.exp_x_interp_calibrated = None
        self.exp_y_interp_calibrated = None
        self.baseline = baseline
        self.calibrated_residual = None
        self.params_to_calibrate = None
        self.raw_layer = raw_layer

    def norm_to(self, file):
        if file is not None:
            self.experiment.norm_to(file=file)

    def slice(self, slice_start=None, slice_end=None):
        self.experiment.slice(slice_start=slice_start, slice_end=slice_end)

    def calibrate(self, source_to_detector_m, offset_us, vary='all', each_step=False):
        """
        calibrate the instrumental parameters: source-to-detector-distance & detector delay
        :param each_step: boolean. True -> show values and chi^2 of each step
        :param source_to_detector_m: estimated distance in m
        :param offset_us: estimated time offset in us
        :param vary: vary one of or both of 'source_to_detector' and 'offset' to calibrate (default: 'all')

        :return: lmfit MinimizerResult
        """
        self.init_source_to_detector_m = source_to_detector_m
        self.init_offset_us = offset_us
        if vary not in ['source_to_detector', 'offset', 'all', 'none']:
            raise ValueError("'vary=' can only be one of ['source_to_detector', 'offset', 'all' 'none']")
        simu_x = self.simu_x
        simu_y = self.simu_y

        source_to_detector_vary_tag = True
        offset_vary_tag = True
        if vary == 'source_to_detector':
            offset_vary_tag = False
        if vary == 'offset':
            source_to_detector_vary_tag = False
        if vary == 'none':
            source_to_detector_vary_tag = False
            offset_vary_tag = False
        self.params_to_calibrate = Parameters()
        self.params_to_calibrate.add('source_to_detector_m', value=source_to_detector_m,
                                     vary=source_to_detector_vary_tag)
        self.params_to_calibrate.add('offset_us', value=offset_us, vary=offset_vary_tag)
        # Print before
        print("Params before calibration:")
        self.params_to_calibrate.pretty_print()
        # Use lmfit to obtain 'source_to_detector_m' & 'offset_us' to minimize 'y_gap_for_calibration'
        self.calibrate_result = minimize(y_gap_for_calibration, self.params_to_calibrate, method='leastsq',
                                         args=(simu_x, simu_y,
                                               self.energy_min, self.energy_max, self.energy_step,
                                               self.experiment, self.baseline, each_step))
        # Print after
        print("Params after calibration:")
        self.calibrate_result.__dict__['params'].pretty_print()
        # Print chi^2
        self.calibrated_residual = self.calibrate_result.__dict__['residual']
        print("Calibration chi^2 : {}\n".format(sum(self.calibrated_residual ** 2)))
        self.calibrated_offset_us = self.calibrate_result.__dict__['params'].valuesdict()['offset_us']
        self.calibrated_source_to_detector_m = \
            self.calibrate_result.__dict__['params'].valuesdict()['source_to_detector_m']

        # Save the calibrated experimental x & y in Calibration class
        self.exp_x_raw_calibrated = self.experiment.x_raw(angstrom=False,
                                                          offset_us=self.calibrated_offset_us,
                                                          source_to_detector_m=self.calibrated_source_to_detector_m)
        self.exp_y_raw_calibrated = self.experiment.y_raw(transmission=False, baseline=self.baseline)

        self.exp_x_interp_calibrated, self.exp_y_interp_calibrated = self.experiment.xy_scaled(
            energy_min=self.energy_min,
            energy_max=self.energy_max,
            energy_step=self.energy_step,
            offset_us=self.calibrated_offset_us,
            source_to_detector_m=self.calibrated_source_to_detector_m,
            baseline=self.baseline)

        return self.calibrate_result

    def plot(self, interp=False, before=False, table=True, grid=False):
        """
        Plot the raw experimental data and theoretical resonance signal after calibration
        :param table: boolean. True -> display table of calibrated parameters
        :param before: boolean. True -> plot the raw data before calibration applied.
        :param interp: boolean. True -> display interpolated exp data
                                False -> display raw exp data
        :return: plot of raw experimental data and theoretical resonance signal after calibration
        """
        simu_label = 'Ideal'
        exp_label = 'Exp'
        exp_before_label = 'Exp_raw'
        exp_interp_label = 'Exp_interp'
        sample_name = ' & '.join(self.layer_list)
        fig_title = 'Calibration result of sample ' + '(' + sample_name + ')'

        # Plot graph
        if table is True:
            ax1 = plt.subplot2grid(shape=(10, 10), loc=(0, 1), rowspan=8, colspan=8)
        else:
            ax1 = plt.subplot(111)
        ax1.plot(self.simu_x, self.simu_y, 'b-', label=simu_label, linewidth=1)
        if before is True:
            ax1.plot(self.experiment.x_raw(offset_us=self.init_offset_us,
                                           source_to_detector_m=self.init_source_to_detector_m),
                     self.experiment.y_raw(baseline=self.baseline),
                     'cs', label=exp_before_label, markersize=2)
        if interp is False:
            ax1.plot(self.exp_x_raw_calibrated, self.exp_y_raw_calibrated, 'rx', label=exp_label, markersize=2)
        else:
            ax1.plot(self.exp_x_interp_calibrated, self.exp_y_interp_calibrated, 'r-.', label=exp_interp_label,
                     linewidth=1)
        ax1.set_xlim([0, self.energy_max])
        ax1.set_ylim(ymax=1.01)
        ax1.set_title(fig_title)
        ax1.set_xlabel('Energy (eV)')
        ax1.set_ylabel('Attenuation')
        ax1.legend(loc='best')
        # ax1.legend(bbox_to_anchor=(1., 1), loc=2, borderaxespad=0.)
        # ax1.legend(bbox_to_anchor=(0, 0.93, 1., .102), loc=3, borderaxespad=0.)
        if grid is True:
            # ax1.set_xticks(np.arange(0, 100, 10))
            # ax1.set_yticks(np.arange(0, 1., 0.1))
            ax1.grid()

        # Plot table
        if table is True:
            # ax2 = plt.subplot2grid(shape=(10, 7), loc=(0, 1), rowspan=4, colspan=5)
            # ax2.axis('off')
            columns = self.calibrate_result.__dict__['var_names']
            rows = ['Before', 'After']
            _row_before = []
            _row_after = []
            for _each in columns:
                _row_after.append(self.calibrate_result.__dict__['params'].valuesdict()[_each])
                _row_before.append(self.params_to_calibrate.valuesdict()[_each])
            table = ax1.table(rowLabels=rows, colLabels=columns,  # colWidths=
                              cellText=[[self.init_source_to_detector_m, self.init_offset_us],
                                        [self.calibrated_source_to_detector_m,
                                         self.calibrated_offset_us]],  # rows of data values
                              bbox=[0, -0.33, 1.0, 0.18]  # [left,bottom,width,height]
                              )
            # table.scale(0.5, 1)
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            plt.tight_layout()

        plt.show()
        # plt.savefig('test.tiff', dpi=300)
