# -*- coding: utf-8 -*-
"""Current defines a generic flux time series that can be associated to a transport coefficient."""

import numpy as np
from .. import md
from ..md.mdsample import MDSample
from ..md.tools.spectrum import freq_THz_to_red, freq_red_to_THz
from ..md.tools.resample import filter_and_sample

from thermocepstrum.utils.loadAfterPlt import plt
from thermocepstrum.utils.utils import PrintMethod
log = PrintMethod()

try:
    plt
except:
    log.write_log('Warning: plt undefined')


class Current(MDSample):
    """
    Current API for thermo-cepstral analysis.
    Defines a Current object with useful tools to perform analysis.

    INPUT:
     - traj          the heat current time series (N * N_COMPONENTS array)
       For a multi-component fluid use a (N_FLUID_COMPONENTS * N * N_COMPONENTS array)
     - UNITS         the units of current ('metal', 'real')
     - DT_FS         MD time step [fs]
     - TEMPERATURE   average temperature [K]
     - VOLUME        simulation cell volume [A^3]
     - PSD_FILTER_W  PSD filter window [freq_units] (optional)
     - FREQ_UNITS    frequency units   [THz or red] (optional)
    """
    _current_type = None
    _input_parameters = {'DT_FS'}
    _optional_parameters = {'PSD_FILTER_W', 'FREQ_UNITS'}

    # parameters are class-specific (a HeatCurrent may use different ones wrt ElectricCurrent)

    def __init__(self, traj, plotter=None, **params):
        # e.g. params: (DT_FS, UNITS, TEMPERATURE, VOLUME, PSD_FILTER_W=None, FREQ_UNITS='THz')
        # validate input parameters
        params = {k.upper(): v for k, v in params.items()}   # convert keys to uppercase
        keyset = set(params.keys())
        if not self._input_parameters.issubset(keyset):
            raise ValueError('The input parameters {} must be defined.'.format(self._input_parameters - keyset))
        if not keyset.issubset(self._input_parameters | self._optional_parameters):
            raise ValueError(
                'The input parameters {} are not valid.'.format(keyset -
                                                                (self._input_parameters | self._optional_parameters)))

        if plotter:
            self.plotter = plotter
        else:
            log.write_log('Warning: plotter not initialized. Plotts will be not created.')
            self.plotter = None

        print(params)

        # pop non unit-specific parameters
        PSD_FILTER_W = params.pop('PSD_FILTER_W', None)
        FREQ_UNITS = params.pop('FREQ_UNITS', 'THz')

        DT_FS = params.get('DT_FS')
        self.initialize_currents(traj, DT_FS)
        self.initialize_units(**params)   # (UNITS, TEMPERATURE, VOLUME, DT_FS)
        if self.traj is not None:
            self.compute_psd(PSD_FILTER_W, FREQ_UNITS)
            self.initialize_cepstral_parameters()
        else:
            log.write_log('Warning: trajectory not initialized. You should manually initialize what you need.')

        self.dct = None

    def __repr__(self):
        msg = type(self).__name__ +\
              '\n  N_CURRENTS  =  {}\n'.format(self.N_CURRENTS)
        for key in self._input_parameters - {'DT_FS'}:
            msg += '  {:11} =  {}\n'.format(key, getattr(self, key))
        msg += super().__repr__()
        if self.otherMD:
            for current in self.otherMD:
                msg += current.__repr__()
        if self.dct:
            msg += self.dct.__repr__()
        return msg

    def _get_builder(self):
        """
        Get a tuple (class, builder) that can be used to build a new object with same parameters:
          TimeSeries, builder = self._get_builder()
          new_ts = TimeSeries(**builder)
        """
        # this is a virtual method
        pass

    def initialize_currents(self, j, DT_FS):
        # check if we have a multicomponent fluid
        j = np.array(j, dtype=float)
        if (len(j.shape) == 3):
            self.N_CURRENTS = j.shape[0]
            if (self.N_CURRENTS == 1):
                self.MANY_CURRENTS = False
                j = np.squeeze(j, axis=0)
            else:
                self.MANY_CURRENTS = True
        elif (len(j.shape) <= 2):
            self.N_CURRENTS = 1
            self.MANY_CURRENTS = False
        else:
            raise ValueError('Shape of j {} not valid.'.format(j.shape))

        if self.MANY_CURRENTS:
            log.write_log('Using multicomponent code.')
            super().__init__(traj=j[0], DT_FS=DT_FS)
            # initialize other MDSample currents
            self.otherMD = [MDSample(traj=js, DT_FS=DT_FS) for js in j[1:]]
        else:
            log.write_log('Using single component code.')
            super().__init__(traj=j, DT_FS=DT_FS)
            self.otherMD = None

    @staticmethod
    def get_units_list():
        # this is a virtual method
        # TODO: find a way to read units from the functions defined in the module 'current/units/*_current_type*.py'
        # TODO: another method should return the function directly from the key
        pass

    def initialize_units(self, **parameters):
        """
        Initializes the units and define the kappa_scale.
        """
        # this is a virtual method
        pass

    def compute_psd(self, PSD_FILTER_W=None, freq_units='THz'):
        # overrides MDSample method
        """
        Compute the periodogram from the heat current time series.
        If a PSD_FILTER_W (expressed in freq_units) is known or given, the psd is also filtered.
        The PSD is multiplied by DT_FS at the end.
        """
        if self.MANY_CURRENTS:
            if self.otherMD is None:
                raise RuntimeError('self.otherMD cannot be None (wrong/missing initialization?)')
            self.compute_kappa_multi(self.otherMD, PSD_FILTER_W, freq_units)
        else:
            super().compute_psd(PSD_FILTER_W, freq_units)

    def initialize_cepstral_parameters(self):
        """
        Defines the parameters of the theoretical distribution of the cepstrum.
        """
        if not self.MANY_CURRENTS:
            self.ck_THEORY_var, self.psd_THEORY_mean = \
                md.cepstral.multicomp_cepstral_parameters(self.NFREQS, self.N_COMPONENTS)
        else:
            if self.ndf_chi is None:
                raise RuntimeError('self.ndf_chi cannot be None.')
            self.ck_THEORY_var, self.psd_THEORY_mean = \
                md.cepstral.multicomp_cepstral_parameters(self.NFREQS, self.ndf_chi)

    def cepstral_analysis(self, aic_type='aic', Kmin_corrfactor=1.0, K_PSD=None):
        """
        Performs Cepstral Analysis on the heat current trajectory.
           aic_type      = the Akaike Information Criterion function used to choose the cutoff ('aic', 'aicc')
           Kmin_corrfactor = correction factor multiplied by the AIC cutoff (cutoff = Kmin_corrfactor * aic_Kmin)

        Resulting conductivity:
            appa_Kmin  +/-  kappa_Kmin_std   [W/(m*K)]
        """

        self.dct = md.CosFilter(self.logpsd, ck_theory_var=self.ck_THEORY_var, \
            psd_theory_mean=self.psd_THEORY_mean, aic_type=aic_type, Kmin_corrfactor=Kmin_corrfactor)
        self.dct.scan_filter_tau(K_PSD=K_PSD)
        self.kappa_Kmin = self.dct.tau_Kmin * self.kappa_scale * 0.5
        self.kappa_Kmin_std = self.dct.tau_std_Kmin * self.kappa_scale * 0.5

        self.cepstral_log = \
              '-----------------------------------------------------\n' +\
              '  CEPSTRAL ANALYSIS\n' +\
              '-----------------------------------------------------\n' +\
              '  AIC_Kmin  = {:d}  (P* = {:d}, corr_factor = {:4f})\n'.format(self.dct.aic_Kmin, self.dct.aic_Kmin + 1, self.dct.Kmin_corrfactor) +\
              '  L_0*   = {:18f} +/- {:10f}\n'.format(self.dct.logtau_Kmin, self.dct.logtau_std_Kmin) +\
              '  S_0*   = {:18f} +/- {:10f}\n'.format(self.dct.tau_Kmin, self.dct.tau_std_Kmin) +\
              '-----------------------------------------------------\n' +\
              '  kappa* = {:18f} +/- {:10f}  W/mK\n'.format(self.kappa_Kmin, self.kappa_Kmin_std) +\
              '-----------------------------------------------------\n'
        log.write_log(self.cepstral_log)

    ###################################
    ###  PLOT METHODS
    ###################################

    def plot_periodogram(self, PSD_FILTER_W=None, freq_units='THz', freq_scale=1.0, axes=None, kappa_units=False,
                         FIGSIZE=None, **plot_kwargs):   # yapf: disable
        if self.check_plotter():
            return self.plotter.plot_periodogram(self, PSD_FILTER_W, freq_units,
                                                 freq_scale, axes, kappa_units,
                                                 FIGSIZE, **plot_kwargs)

    def plot_ck(self, axes=None, label=None, FIGSIZE=None):
        if self.check_plotter():
            return self.plotter.plot_ck(self, axes, label, FIGSIZE)

    def plot_L0_Pstar(self, axes=None, label=None, FIGSIZE=None):
        if self.check_plotter():
            return self.plot_L0_Pstar(self, axes, label, FIGSIZE)

    def plot_kappa_Pstar(self, axes=None, label=None, FIGSIZE=None):
        if self.check_plotter():
            return self.plotter.plot_kappa_Pstar(self, axes, label, FIGSIZE)

    def plot_cepstral_spectrum(self, freq_units='THz', freq_scale=1.0, axes=None, kappa_units=True, FIGSIZE=None,
                               **plot_kwargs):   # yapf: disable
        if self.check_plotter():
            return self.plot_cepstral_spectrum(self, freq_units, freq_scale, axes, kappa_units, FIGSIZE, **plot_kwargs)

    def resample(self, TSKIP=None, fstar_THz=None, FILTER_W=None, plot=True, PSD_FILTER_W=None,
                 freq_units='THz', FIGSIZE=None, verbose=True):   # yapf: disable
        """
        Simulate the resampling of the time series.

        Parameters
        ----------
        TSKIP        = sampling time [steps]
        fstar_THz    = target cutoff frequency [THz]
        TSKIP and fstar_THZ are mutually exclusive.

        FILTER_W     = pre-sampling filter window width [steps]
        plot         = plot the PSD [True]
        PSD_FILTER_W = PSD filtering window width [chosen frequency units]
        freq_units   = 'thz'  [THz]
                       'red'  [omega*DT/(2*pi)]
        FIGSIZE      = plot figure size
        verbose      = print log [True]

        Returns
        -------
        xf : a filtered & resampled time series object
        ax : an array of plot axes, optional (if plot=True)
        """
        xf = super().resample(TSKIP, fstar_THz, FILTER_W, False, PSD_FILTER_W, freq_units, None, verbose)

        if plot and self.check_plotter():
            return self.plotter.plt_resample(self, xf, freq_units, PSD_FILTER_W, FIGSIZE)
        return xf

    def fstar_analysis(self, TSKIP_LIST, aic_type='aic', Kmin_corrfactor=1.0, plot=True, axes=None, FIGSIZE=None,
                       verbose=False, **plot_kwargs):   # yapf: disable
        return fstar_analysis(self, TSKIP_LIST, aic_type, Kmin_corrfactor, plot, axes, FIGSIZE, verbose, **plot_kwargs)

    def check_plotter(self):
        if self.plotter:
            return True
        else:
            return False


################################################################################


def fstar_analysis(x, TSKIP_LIST, aic_type='aic', Kmin_corrfactor=1.0, plot=True, axes=None, FIGSIZE=None,
                   verbose=False, **plot_kwargs):   # yapf: disable
    """
    Perform cepstral analysis on a set of resampled time series, to study the effect of f*.
    For each TSKIP in TSKIP_LIST, the HeatCurrent x is filtered & resampled, and then cesptral-analysed.

    Parameters
    ----------
    TSKIP_LIST    = list of sampling times [steps]
    aic_type      = the Akaike Information Criterion function used to choose the cutoff ('aic', 'aicc')
    Kmin_corrfactor = correction factor multiplied by the AIC cutoff (cutoff = Kmin_corrfactor * aic_Kmin)
    plot          = plot the PSD (True/False)
    axes          = matplotlib.axes.Axes object (if None, create one)
    FIGSIZE       = plot figure size
    verbose       = verbose output (default: False)
    **plot_kwargs = other parameters passed to plot function

    Returns
    -------
    xf : array_like
        an array of HeatCurrents, corresponding the values of TSKIP_LIST
    ax : array_like, optional (if plot=True)
        array of plot axes
    fig : figure, optional (if axes=None)
        plot figure
    """

    if not isinstance(x, HeatCurrent):
        raise ValueError('x must be a HeatCurrent object.')

    xf = []
    for TSKIP in TSKIP_LIST:
        log.write_log('TSKIP =  {:d}'.format(TSKIP))
        xff = resample_current(x, TSKIP, plot=False, verbose=verbose)
        xff.cepstral_analysis(aic_type, Kmin_corrfactor)
        xf.append(xff)
    FSTAR_THZ_LIST = [xff.Nyquist_f_THz for xff in xf]

    if plot and x.check_plotter():
        return x.plotter.plot_fstar_analysis(x, xf, FSTAR_THZ_LIST, axes, FIGSIZE, **plot_kwargs)
    else:
        return xf
