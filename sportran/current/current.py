# -*- coding: utf-8 -*-
"""Current defines a generic flux time series that can be associated to a transport coefficient."""

import numpy as np
import inspect
from sportran import md
from sportran.md.mdsample import MDSample
from sportran.md.tools.spectrum import freq_THz_to_red, freq_red_to_THz
from sportran.md.tools.resample import filter_and_sample
from . import units
from sportran.utils import log
from sportran.plotter.current import CurrentPlotter

__all__ = ('Current', 'fstar_analysis',)


class Current(MDSample):
    """
    Current API for thermo-cepstral analysis.
    Defines a Current object with useful tools to perform analysis.

    INPUT parameters:
     - traj          the heat current time series array (N * N_COMPONENTS array)
       For a multi-component fluid use a (N_FLUID_COMPONENTS * N * N_COMPONENTS array)
     - DT_FS         MD time step [fs]
     - KAPPA_SCALE   the GK conversion factor, multiplies by the GK integral

    OPTIONAL parameters:
     - PSD_FILTER_W  PSD filter window [freq_units] (optional)
     - FREQ_UNITS    frequency units   [THz or red] (optional)
     - MAIN_CURRENT_INDEX for a multi-current time series, the index of the "main" current (e.g. energy) [0]
     - MAIN_CURRENT_FACTOR factor to be multiplied by the main current [1.0]

    The default plotter is `plotter.CurrentPlotter`. It can be set by `Current.set_plotter`.
    """

    # parameters are class-specific (a HeatCurrent may use different ones wrt ElectricCurrent) and case-insensitive
    _current_type = None
    _input_parameters = {'DT_FS', 'KAPPA_SCALE'}
    _optional_parameters = {'PSD_FILTER_W', 'FREQ_UNITS', 'MAIN_CURRENT_INDEX', 'MAIN_CURRENT_FACTOR'}
    _KAPPA_SI_UNITS = ''
    _default_plotter = CurrentPlotter

    def __init__(self, traj, **params):
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

        # pop non unit-specific parameters
        PSD_FILTER_W = params.pop('PSD_FILTER_W', None)
        FREQ_UNITS = params.pop('FREQ_UNITS', 'THz')

        DT_FS = params.pop('DT_FS')
        MAIN_CURRENT_INDEX = params.pop('MAIN_CURRENT_INDEX', 0)
        MAIN_CURRENT_FACTOR = params.pop('MAIN_CURRENT_FACTOR', 1.0)
        self.initialize_currents(traj, DT_FS, MAIN_CURRENT_INDEX, MAIN_CURRENT_FACTOR)
        self.initialize_units(**params)   # KAPPA_SCALE or (e.g. UNITS, TEMPERATURE, VOLUME)
        if self.traj is not None:
            self.compute_psd(PSD_FILTER_W, FREQ_UNITS)
            self.initialize_cepstral_parameters()
        else:
            log.write_log('Warning: trajectory not initialized. You should manually initialize what you need.')
        self.dct = None

    def __repr__(self):
        msg = type(self).__name__ +\
              '\n  N_CURRENTS  =  {}\n'.format(self.N_CURRENTS) +\
              '  KAPPA_SCALE =  {}\n'.format(self.KAPPA_SCALE)
        for key in self._input_parameters - {'DT_FS', 'KAPPA_SCALE'}:
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
        # TODO: unify this method across the subclasse
        if self.MANY_CURRENTS:
            traj_array = np.row_stack(([self.traj], [j.traj for j in self.otherMD]))
        else:
            traj_array = self.traj
        builder = dict(traj=traj_array, DT_FS=self.DT_FS, KAPPA_SCALE=self.KAPPA_SCALE,
                       PSD_FILTER_W=self.PSD_FILTER_W_THZ, FREQ_UNITS='THz')
        return type(self), builder

    @classmethod
    def set_plotter(cls, plotter=None):
        """
        Set the plotter class.
        The _plotter attribute will contain the selected plotter class.
        All the plot functions of plotter (named 'plot_*') will be transformed into methods of Current.

        **NOTE**
        If called by a subclass, it will change the plotter of the base class (Current) and all its subclasses.
        If this is not a good behavior, we should change it in the future.
        """
        # if called by a subclass of Current, change the base class (Current)
        if issubclass(cls, Current) and cls != Current:
            cls = Current
        # (note: it is not possible to delete the parent class' attributes from a child. But here we forcibly do this operation on cls = Current)
        super().set_plotter(plotter)

    def initialize_currents(self, j, DT_FS, main_current_index=0, main_current_factor=1.0):
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
            super().__init__(traj=(j[main_current_index] * main_current_factor), DT_FS=DT_FS)
            # initialize other MDSample currents
            other_currents_idxs = (np.arange(self.N_CURRENTS) != main_current_index)   # select the other currents
            self.otherMD = [MDSample(traj=js, DT_FS=DT_FS) for js in j[other_currents_idxs]]
        else:
            log.write_log('Using single component code.')
            super().__init__(traj=(j * main_current_factor), DT_FS=DT_FS)
            self.otherMD = None

    @classmethod
    def _get_units(cls):
        try:
            # get the units submodule corresponding to this class
            units_module = getattr(units, cls._current_type)
        except AttributeError:
            print('No units submodule defined for the current type "{}". Add units to a file "current/units/{}.py".'.
                  format(cls._current_type, cls._current_type))
            return {}
        except TypeError:
            raise RuntimeError('No units can be defined for a generic Current. Define a "KAPPA_SCALE" instead.')

        # get all functions that start with "scale_kappa_" into a dictionary {"name": function}
        units_prefix = 'scale_kappa_'
        units_d = {
            name.replace(units_prefix, ''): function for name, function in inspect.getmembers(
                units_module, predicate=lambda f: inspect.isfunction(f) and f.__name__.startswith(units_prefix))
        }
        if not units_d:
            print(
                'Warning: No units defined for a current type "{}". Add them to the module "current/units/{}.py'.format(
                    cls._current_type, cls._current_type))
        return units_d

    @classmethod
    def get_units_list(cls):
        """
        Get the list of supported units.
        Units are defined in the module current/units/{current_type}.py, where
        {current_type} is the _current_type attribute of this class ('heat', 'electric', ...).
        """
        return cls._get_units().keys()

    def initialize_units(self, **parameters):
        """
        Initializes the units and defines the KAPPA_SCALE.
        """
        self.UNITS = parameters.pop('UNITS', None)

        # set unit-specific parameters
        for param, value in parameters.items():
            self.__setattr__(param, value)

        # validate units and define KAPPA_SCALE from units conversion function
        if self.UNITS:
            units_list = self.get_units_list()
            if len(units_list) == 0:
                raise RuntimeError(
                    'No units defined for a current type "{}". Add them to the module "current/units/{}.py'.format(
                        self._current_type, self._current_type))
            elif self.UNITS in units_list:
                units_conversion_func = self._get_units()[self.UNITS]
                self.KAPPA_SCALE = units_conversion_func(**parameters)
            else:
                raise ValueError('Units "{}" not valid. Valid units are:\n  {}'.format(
                    self.UNITS, self.get_units_list()))

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
           kappa_Kmin  +/-  kappa_Kmin_std   [SI units]
        """

        self.dct = md.CosFilter(self.logpsd, ck_theory_var=self.ck_THEORY_var, \
            psd_theory_mean=self.psd_THEORY_mean, aic_type=aic_type, Kmin_corrfactor=Kmin_corrfactor)
        self.dct.scan_filter_tau(K_PSD=K_PSD)
        self.kappa_Kmin = self.dct.tau_Kmin * self.KAPPA_SCALE * 0.5
        self.kappa_Kmin_std = self.dct.tau_std_Kmin * self.KAPPA_SCALE * 0.5

        self.cepstral_log = \
              '-----------------------------------------------------\n' +\
              '  CEPSTRAL ANALYSIS\n' +\
              '-----------------------------------------------------\n' +\
              '  AIC_Kmin  = {:d}  (P* = {:d}, corr_factor = {:4f})\n'.format(self.dct.aic_Kmin, self.dct.aic_Kmin + 1, self.dct.Kmin_corrfactor) +\
              '  L_0*   = {:18f} +/- {:10f}\n'.format(self.dct.logtau_Kmin, self.dct.logtau_std_Kmin) +\
              '  S_0*   = {:18f} +/- {:10f}\n'.format(self.dct.tau_Kmin, self.dct.tau_std_Kmin) +\
              '-----------------------------------------------------\n' +\
              '  kappa* = {:18f} +/- {:10f}  {}\n'.format(self.kappa_Kmin, self.kappa_Kmin_std, self._KAPPA_SI_UNITS) +\
              '-----------------------------------------------------\n'
        log.write_log(self.cepstral_log)

    def resample(self, TSKIP=None, fstar_THz=None, FILTER_W=None, plot=False, PSD_FILTER_W=None,
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

        if plot:
            try:
                axs = self.plot_resample(xf=xf, freq_units=freq_units, PSD_FILTER_W=PSD_FILTER_W, FIGSIZE=FIGSIZE)
                return xf, axs
            except AttributeError:
                print('Plotter does not support the plot_resample method')
        else:
            return xf

    def fstar_analysis(self, TSKIP_LIST, aic_type='aic', Kmin_corrfactor=1.0, plot=True, axes=None, FIGSIZE=None,
                       verbose=False, **plot_kwargs):   # yapf: disable
        return fstar_analysis(self, TSKIP_LIST, aic_type, Kmin_corrfactor, plot, axes, FIGSIZE, verbose, **plot_kwargs)


################################################################################

# set the default plotter of this class
Current.set_plotter()

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

    if not isinstance(x, Current):
        raise ValueError('x must be a Current object or a subclass.')

    xf = []
    for TSKIP in TSKIP_LIST:
        log.write_log('TSKIP = {:4d} - FSTAR = {:8g} THz'.format(TSKIP, x.Nyquist_f_THz / TSKIP))
        xff = x.resample(TSKIP=TSKIP, plot=False, verbose=verbose)
        xff.cepstral_analysis(aic_type, Kmin_corrfactor)
        xf.append(xff)
    FSTAR_THZ_LIST = [xff.Nyquist_f_THz for xff in xf]

    if plot:
        try:
            x.plot_fstar_analysis(xf, FSTAR_THZ_LIST, axes=axes, FIGSIZE=FIGSIZE, **plot_kwargs)
        except AttributeError:
            print('Plotter does not support the plot_resample method')
    else:
        return xf
