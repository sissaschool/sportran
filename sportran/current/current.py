# -*- coding: utf-8 -*-
"""Current defines a generic flux time series that can be associated to a transport coefficient."""

import numpy as np
import inspect
from sportran.md.mdsample import MDSample
from sportran.md.cepstral import CepstralFilter, multicomp_cepstral_parameters
from sportran.md.tools.filter import runavefilter
from sportran.md.tools.spectrum import freq_THz_to_red, freq_red_to_THz
from sportran.md.tools.resample import filter_and_sample
from . import units
from sportran.utils import log
from sportran.plotter.current import CurrentPlotter

__all__ = ['Current', 'fstar_analysis']


class Current(MDSample):
    """
    Current API for thermo-cepstral analysis.
    Defines a Current object with useful tools to perform analysis.

    INPUT parameters:
     - traj          the current time series (N, N_EQUIV_COMPONENTS) array
       For a multi-component fluid use a (N_CURRENTS, N, N_EQUIV_COMPONENTS) array
     - DT_FS         MD time step [fs]
     - KAPPA_SCALE   the GK conversion factor, multiplies the GK integral

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
        self.cepf = None

    def __repr__(self):
        msg = type(self).__name__ +\
              '\n  N_CURRENTS  =  {}\n'.format(self.N_CURRENTS) +\
              '  KAPPA_SCALE =  {}\n'.format(self.KAPPA_SCALE)
        for key in self._input_parameters - {'DT_FS', 'KAPPA_SCALE'}:
            msg += '  {:11} =  {}\n'.format(key, getattr(self, key))
        msg += super().__repr__()
        if self.otherMD:
            msg += 'additional currents:\n'
            for current in self.otherMD:
                msg += current.__repr__()
        if self.cepf:
            msg += self.cepf.__repr__()
        try:
            msg += '\n  kappa* = {:18f} +/- {:10f}  {}\n'.format(self.kappa, self.kappa_std, self._KAPPA_SI_UNITS)
        except AttributeError:
            pass
        return msg

    @property
    def _builder(self):
        """
        Returns a dictionary of all keyworded parameters needed to rebuild an identical object of the same class.
        The trajectory is excluded. Used by self._get_builder().
        """
        return dict(DT_FS=self.DT_FS, KAPPA_SCALE=self.KAPPA_SCALE, PSD_FILTER_W=self.PSD_FILTER_W_THZ,
                    FREQ_UNITS='THz')

    def _get_builder(self):
        """
        Get a tuple (class, builder) that can be used to build a new object with same parameters:
          TimeSeries, builder = self._get_builder()
          new_ts = TimeSeries(**builder)
        """
        if self.MANY_CURRENTS:
            traj_array = np.row_stack(([self.traj], [j.traj for j in self.otherMD]))
        else:
            traj_array = self.traj
        kwargs = self._builder
        kwargs.update(traj=traj_array)
        return type(self), kwargs

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

        # initialize cospectrum, that is not initialized by MDSample.initialize_psd
        self.cospectrum = None
        self.fcospectrum = None

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
        # number of degrees of freedom of the chi-square distribution of the psd / 2
        self.ndf_chi = self.N_EQUIV_COMPONENTS - self.N_CURRENTS + 1
        if self.ndf_chi <= 0:
            raise RuntimeError('The number of degrees of freedom of the chi-squared distribution is <=0. The number of '
                               'equivalent (Cartesian) components of the input current must be >= number of currents.')

        if self.MANY_CURRENTS:
            if self.otherMD is None:
                raise RuntimeError('self.otherMD cannot be None (wrong/missing initialization?)')
            self._compute_psd_multi(self.otherMD, PSD_FILTER_W, freq_units)
        else:
            super().compute_psd(PSD_FILTER_W, freq_units)

    def _compute_psd_multi(self, others, PSD_FILTER_W=None, freq_units='THz', normalize=False, call_other=True):
        """
        For multi-component (many-current) systems: compute the cospectrum matrix and the transport coefficient.
        The results have almost the same statistical properties.
        The chi-square distribution has ndf = 2(ndf_chi) = 2(l - M + 1), where l is the number of time series for each
        current (N_EQUIV_COMPONENTS), M is the number of currents (N_CURRENTS).
        ! NOTICE: if l < M this will not work.

        In this routine the mean over the number of temporal series is already multiplied by the correct factor (the
        transport coefficient will be obtained by multiplying the result by 0.5, as in the one-component case).
        The output arrays are the same as in the one-component case.
        The elements of the matrix are multiplied by DT_FS at the end.
        If a PSD_FILTER_W is known or given, the psd is also filtered.
        others is a list of other currents, i.e. MDSample objects.

        For example, in the case of 4 currents, of which j is the energy current and j1, j2, j3 are mass currents:
           j._compute_psd_multi([j1,j2,j3], PSD_FILTER_W, freq_units)
        """
        # check if others is an array
        if not isinstance(others, (list, tuple, np.ndarray)):
            others = [others]
        if self.traj is None:
            raise ValueError('Trajectory not defined.')

        self.spectrALL = np.fft.rfft(self.traj, axis=0)
        self.NFREQS = self.spectrALL.shape[0]
        self.freqs = np.linspace(0., 0.5, self.NFREQS)
        self.DF = 0.5 / (self.NFREQS - 1)
        self.DF_THZ = freq_red_to_THz(self.DF, self.DT_FS)
        self.freqs_THz = self.freqs / self.DT_FS * 1000.
        self.Nyquist_f_THz = self.freqs_THz[-1]

        # calculate the same thing on the other trajectory
        if (call_other):
            for other in others:   # call other._compute_psd_multi (MDsample method)
                Current._compute_psd_multi(other, [self], PSD_FILTER_W, freq_units, normalize, False)
        else:
            return

        # define the cospectrum matrix. Its shape is (2, 2, NFREQS, n_spatial_dim)
        #  [  self.spectrALL*self.spectrALL.conj()     self.spectrALL*other.spectrALL.conj() ]
        #  [ other.spectrALL*self.spectrALL.conj()    other.spectrALL*other.spectrALL.conj() ]
        other_spectrALL = []
        for other in others:
            other_spectrALL.append(other.spectrALL)

        # compute the matrix defined by the outer product of only the first indexes of the two arrays
        covarALL = self.DT_FS / (2. * (self.NFREQS - 1.)) *\
                    np.einsum('a...,b...->ab...', np.array([self.spectrALL] + other_spectrALL),
                                                  np.array([self.spectrALL] + other_spectrALL).conj())

        # number of degrees of freedom of the chi-square distribution of the psd / 2
        assert self.ndf_chi == (covarALL.shape[3] - len(other_spectrALL))

        # compute the sum over the last axis (equivalent Cartesian components):
        self.cospectrum = covarALL.sum(axis=3)

        # compute the element 1/"(0,0) of the inverse" (aka the transport coefficient)
        # the diagonal elements of the inverse have very convenient statistical properties
        multi_psd = (np.linalg.inv(self.cospectrum.transpose((2, 0, 1)))[:, 0, 0]**-1).real / self.ndf_chi

        if normalize:
            multi_psd = multi_psd / np.trapz(multi_psd) / self.N / self.DT_FS

        self.psd = multi_psd
        self.logpsd = np.log(self.psd)
        self.psd_min = np.min(self.psd)
        self.psd_power = np.trapz(self.psd)   # one-side PSD power
        if (PSD_FILTER_W is not None) or (self.PSD_FILTER_W is not None):
            self.filter_psd(PSD_FILTER_W, freq_units)

    def filter_psd(self, PSD_FILTER_W=None, freq_units='THz', window_type='rectangular', logpsd_filter_type=1):
        """
        Filter the periodogram with the given PSD_FILTER_W [freq_units].
          - PSD_FILTER_W  PSD filter window [freq_units]
          - freq_units    frequency units   ['THz', 'red' (default)]
          - window_type   filtering window type ['rectangular']
        """
        super().filter_psd(PSD_FILTER_W, freq_units, window_type, logpsd_filter_type)

        if (window_type == 'rectangular'):
            # try to filter the other currents (if present)
            if self.cospectrum is not None:
                self.fcospectrum = []
                for i in range(self.cospectrum.shape[0]):
                    self.fcospectrum.append([])
                    for j in range(self.cospectrum.shape[1]):
                        ffpsd = runavefilter(self.cospectrum[i, j], self.PSD_FILTER_WF)
                        self.fcospectrum[i].append(ffpsd / self.N_EQUIV_COMPONENTS)
                self.fcospectrum = np.asarray(self.fcospectrum)

    def initialize_cepstral_parameters(self):
        """
        Defines the parameters of the theoretical distribution of the cepstrum.
        """
        if not self.MANY_CURRENTS:
            self.ck_THEORY_var, self.psd_THEORY_mean = multicomp_cepstral_parameters(
                self.NFREQS, self.N_EQUIV_COMPONENTS)
        else:
            if self.ndf_chi is None:
                raise RuntimeError('self.ndf_chi cannot be None.')
            self.ck_THEORY_var, self.psd_THEORY_mean = multicomp_cepstral_parameters(self.NFREQS, self.ndf_chi)

    def cepstral_analysis(self, aic_type='aic', aic_Kmin_corrfactor=1.0, manual_cutoffK=None):
        """
        Performs Cepstral Analysis on the Current's trajectory.

        `cutoffK` = (P*-1) is the number of cepstral coefficients retained by the filter.
        By default, this is chosen as the number of cepstral coefficients that minimizes the Akaike Information Criterion,
        multiplied by a correction factor (`aic_Kmin_corrfactor`):
           self.cfilt.cutoffK = argmin(self.cfilt.aic) * aic_Kmin_corrfactor
        This choice can be manually overridden by setting `manual_cutoffK` to the desired value.

        Input parameters:
           aic_type            = the Akaike Information Criterion function used to choose the cutoff ('aic', 'aicc')
           aic_Kmin_corrfactor = correction factor multiplied by the AIC cutoff (cutoffK = aic_Kmin * Kmin_corrfactor) (default: 1.0)
           manual_cutoffK      = (P*-1) = manual cutoff. If set, the AIC cutoff will be ignored.

        The resulting conductivity is returned in the chosen units, generally:
            kappa  +/-  kappa_std   [SI units]

        The log of the analysis can be retried from the variable `self.cepstral_log`.
        """

        self.cepf = CepstralFilter(self.logpsd, ck_theory_var=self.ck_THEORY_var, \
            psd_theory_mean=self.psd_THEORY_mean, aic_type=aic_type)
        self.cepf.scan_filter_tau(cutoffK=manual_cutoffK, aic_Kmin_corrfactor=aic_Kmin_corrfactor)
        self.kappa = self.cepf.tau_cutoffK * self.KAPPA_SCALE * 0.5
        self.kappa_std = self.cepf.tau_std_cutoffK * self.KAPPA_SCALE * 0.5

        self.cepstral_log = \
              '-----------------------------------------------------\n' +\
              '  CEPSTRAL ANALYSIS\n' +\
              '-----------------------------------------------------\n'
        if not self.cepf.manual_cutoffK_flag:
            self.cepstral_log += \
                '  cutoffK = (P*-1) = {:d}  (auto, AIC_Kmin = {:d}, corr_factor = {:4})\n'.format(self.cepf.cutoffK, self.cepf.aic_Kmin, self.cepf.aic_Kmin_corrfactor)
        else:
            self.cepstral_log += \
                '  cutoffK  = (P*-1) = {:d}  (manual, AIC_Kmin = {:d})\n'.format(self.cepf.cutoffK, self.cepf.aic_Kmin, self.cepf.aic_Kmin_corrfactor)
        self.cepstral_log += \
              '  L_0*   = {:18f} +/- {:10f}\n'.format(self.cepf.logtau_cutoffK, self.cepf.logtau_std_cutoffK) +\
              '  S_0*   = {:18f} +/- {:10f}\n'.format(self.cepf.tau_cutoffK, self.cepf.tau_std_cutoffK) +\
              '-----------------------------------------------------\n' +\
              '  kappa* = {:18f} +/- {:10f}  {}\n'.format(self.kappa, self.kappa_std, self._KAPPA_SI_UNITS) +\
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

    def fstar_analysis(self, TSKIP_LIST, aic_type='aic', aic_Kmin_corrfactor=1.0, manual_cutoffK=None, plot=True,
                       axes=None, FIGSIZE=None, verbose=False, **plot_kwargs):   # yapf: disable
        return fstar_analysis(self, TSKIP_LIST, aic_type, aic_Kmin_corrfactor, manual_cutoffK, plot, axes, FIGSIZE,
                              verbose, **plot_kwargs)


################################################################################

# set the default plotter of this class
Current.set_plotter()

################################################################################


def fstar_analysis(x, TSKIP_LIST, aic_type='aic', aic_Kmin_corrfactor=1.0, manual_cutoffK=None, plot=True, axes=None,
                   FIGSIZE=None, verbose=False, **plot_kwargs):   # yapf: disable
    """
    Perform cepstral analysis on a set of resampled time series, to study the effect of f*.
    For each TSKIP in TSKIP_LIST, the HeatCurrent x is filtered & resampled, and then cesptral-analysed.

    Parameters
    ----------
    TSKIP_LIST    = list of sampling times [steps]
    aic_type      = the Akaike Information Criterion function used to choose the cutoff ('aic', 'aicc')
    aic_Kmin_corrfactor = correction factor multiplied by the AIC cutoff (cutoffK = aic_Kmin * aic_Kmin_corrfactor)
    manual_cutoffK = (P*-1) = manual cutoff. If set, the AIC cutoff will be ignored.

    plot          = plot the PSD (default: True)
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
        xff.cepstral_analysis(aic_type=aic_type, aic_Kmin_corrfactor=aic_Kmin_corrfactor, manual_cutoffK=manual_cutoffK)
        xf.append(xff)
    FSTAR_THZ_LIST = [xff.Nyquist_f_THz for xff in xf]

    if plot:
        from sportran.plotter.plotter import plot_fstar_analysis
        try:
            return plot_fstar_analysis(xf, FSTAR_THZ_LIST, original_current=x, axes=axes, FIGSIZE=FIGSIZE,
                                       **plot_kwargs)
        except AttributeError:
            print('Plotter does not support the plot_resample method')
    else:
        return xf
