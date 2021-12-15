# -*- coding: utf-8 -*-

import numpy as np
from scipy.signal import periodogram
from .tools.spectrum import freq_THz_to_red, freq_red_to_THz
from .tools.filter import runavefilter
from .tools.acf import acovf, integrate_acf
from .resample import resample_timeseries
from sportran.utils import log
from sportran.utils.decorators import add_method
from sportran.plotter import Plotter, use_plot_style
from sportran.plotter.mdsample import MDSamplePlotter

__all__ = ('MDSample',)


class MDSample(object):
    """
    An MDSample object contains all the information that represent a
    single and unique Molecular Dynamics sample.
    For example it may contain:
     - a trajectory (any N-dim time series in real space)
     - its spectrum (the Fourier transform)
     - its Power Spectral Density (aka periodogram)
     - ...
    All the information contained in this object is always consistent,
    i.e. it represents a single MD sample. Any operation that alters
    any of the sample's properties should create a new MDSample object, in
    order to preserve the 1:1 correspondence among the sample's attributes.

    An MDSample object can be initialized from any of its main properties
    (trajectory, spectrum, periodogram), although e.g. it is not be possible to
    uniquely define a trajectory from its periodogram, as this contains less
    information.

    ATTRIBUTES:
       - traj       the trajectory (any n-dim real time series).
                    It contains N points.
       - spectr     the spectrum, i.e. trajectory's FFT.
                    For now it is assumed to be one sided, i.e. it contains
                    NFREQS = N/2+1 normalized frequencies contained
                    in the interval [0, 1/(2N*DT)]
       - psd        the Power Spectral Density (periodogram), defined as
                              DT    N-1
                     I(f) =  ---- * SUM | x[n] * exp(-2.0J*pi*f/N) |^2
                               N    n=0
                    with f = [0, 1/(2N*DT)]
       - N          size of traj
       - NFREQS     number of trajectories, should be N/2+1
       - freqs      an array of frequencies, should be [0, 1/(2N*DT)]
       - freqs_THz  an array of frequencies, expressed in THz

    MEMBERS:
        self.DT_FS                  timestep in femtoseconds
        self.fpsd                   filtered periodogram
        self.flogpsd                filtered log-periodogram
        self.acf                    autocorrelation function
        self.N_COMPONENTS           number of EQUIVALENT (cartesian) components (an average over them will be computed)
        self.MANY_COMPONENTS        True if N_COMPONENTS > 1
        self.PSD_FILTER_W           width of the moving average filter (reduced frequency units)
        self.PSD_FILTER_W_THZ       width of the moving average filter (THz)
        self.PSD_FILTER_WF          width of the moving average filter (number of frequencies)

    """

    _default_plotter = MDSamplePlotter

    def __init__(self, traj=None, spectr=None, psd=None, freqs=None, DT_FS=1.0):
        self.DT_FS = DT_FS
        self.initialize_traj(traj)
        self.initialize_spectrum(spectr)
        self.initialize_psd(freqs=freqs, psd=psd, DT_FS=DT_FS)

        self.cospectrum = None
        self.fcospectrum = None

    def __repr__(self):
        msg = 'MDSample:\n' + \
              '  DT_FS:  {}  fs\n'.format(self.DT_FS) + \
              '  traj:   {}  steps  *  {} components\n'.format(self.N, self.N_COMPONENTS) + \
              '          {}  fs\n'.format(None if self.traj is None else self.DT_FS * self.N)
        if self.spectr is not None:
            msg += '  spectr: {}  frequencies\n'.format(self.NFREQS)
        if self.psd is not None:
            msg += '  psd:    {}  frequencies\n'.format(self.psd.size) + \
                   '      DF =      {}  [omega*DT/(2*pi)]\n'.format(self.DF) + \
                   '                {}  [THz]\n'.format(self.DF_THZ) + \
                   '      Nyquist Frequency = {}  [THz]\n'.format(self.Nyquist_f_THz)
        if self.fpsd is not None:
            msg += '  fpsd:   {}  frequencies\n'.format(self.fpsd.size) +\
                   '      PSD_FILTER_W  = {} [omega*DT/(2*pi)]\n'.format(self.PSD_FILTER_W) +\
                   '                    = {} [THz]\n'.format(self.PSD_FILTER_W_THZ) +\
                   '      PSD_FILTER_WF = {} frequencies\n'.format(self.PSD_FILTER_WF)
        if self.acf is not None:
            msg += '  acf:    {}  lags\n'.format(self.NLAGS)
        return msg

    def _get_builder(self):
        """
        Get a tuple (class, builder) that can be used to build a new object with same parameters:
          TimeSeries, builder = self._get_builder()
          new_ts = TimeSeries(**builder)
        """
        builder = dict(traj=self.traj, DT_FS=self.DT_FS)
        return type(self), builder

    @classmethod
    def set_plotter(cls, plotter=None):
        """
        Set the plotter class.
        The _plotter attribute will contain the selected plotter class.
        All the plot functions of plotter (named 'plot_*') will be transformed into methods of Current.
        """
        if plotter is None:
            plotter = cls._default_plotter
        if not (isinstance(plotter, Plotter) or issubclass(plotter, Plotter)):
            raise TypeError('Invalid plotter')

        cls._plotter = plotter
        use_plot_style(plotter._plot_style)

        # delete any plot function already present in this class
        for funcname in filter(lambda name: name.startswith('plot_'),
                               dir(cls)):   # same as [name for name in dir(plotter) if name.startswith('plot_')
            obj = getattr(cls, funcname)
            if callable(obj):
                #print('deleting {} from class {}'.format(obj, cls))
                try:
                    delattr(cls, funcname)
                except AttributeError:
                    pass

        # loop over all functions of the plotter class, and transform them into methods of Current
        for funcname in filter(lambda name: name.startswith('plot_'), dir(plotter)):
            obj = getattr(plotter, funcname)
            if callable(obj):
                add_method(cls)(obj)
                #print('{} added to class {}'.format(obj, cls))

    #############################################
    ###################################
    ###  INITIALIZE METHODS
    ###################################
    #############################################

    def initialize_traj(self, array):
        """
        Initialize a trajectory from an array.
        The dimensions of the array should be:
          (number of time points, number of equivalent components)
        or, in the case of 1 component:
          (number of time points)
        """
        if not isinstance(array, (list, np.ndarray, tuple)):
            raise TypeError('Input trajectory must be an array.')
        if array is not None:
            array = np.array(array, dtype=float)
            if (len(array.shape) == 1):
                self.MANY_COMPONENTS = False
                self.traj = array[:, np.newaxis]
            elif (len(array.shape) == 2):
                self.MANY_COMPONENTS = True
                if (array.shape[0] % 2 == 1):
                    self.traj = array[:-1]
                    log.write_log('Trajectory has an odd number of points. Removing the last one.')
                else:
                    self.traj = array
            else:
                raise TypeError('Input trajectory array has > 2 dimensions.')
            self.N, self.N_COMPONENTS = self.traj.shape
            if (self.N < 2) or (self.N_COMPONENTS < 1):
                raise ValueError('Input trajectory size too small (N = {}, N_COMPONENTS = {}).'.format(
                    self.N, self.N_COMPONENTS))
        else:
            self.traj = None
            self.N = None
            self.N_COMPONENTS = None
        self.acf = None
        self.NLAGS = None

    def initialize_spectrum(self, array):
        if array is not None:
            self.spectr = np.array(array, dtype=complex)
            self.NFREQS = self.spectr.size
            self.DF = 0.5 / (self.NFREQS - 1)
        else:
            self.spectr = None
            self.NFREQS = None
            self.DF = None

    def initialize_psd(self, freq_psd=None, psd=None, freqs=None, DT_FS=None):
        """
        Initialize the PSD. This can be done in 3 ways:
          - passing a tuple  (freqs, psd)
              e.g.   initialize_psd((freqs,psd))
          - passing frequencies and PSD separately
              e.g.   initialize_psd(freqs, psd)
              e.g.   initialize_psd(freqs=freqs, psd=psd)
          - passing PSD only (frequencies will be computed automatically)
              e.g.   initialize_psd(psd)
        """
        # frequencies
        if freq_psd is not None:   # use freq_psd variable
            if (len(freq_psd) == 2):   # (freqs, psd) tuple was passed
                if (freqs is not None) or (psd is not None):
                    raise ValueError('Too many arguments.')
                frequencies = freq_psd[0]
                array = freq_psd[1]
            elif (len(freq_psd) > 2):   # array used as psd or freqs
                if psd is None:   # only psd was passed
                    if freqs is not None:
                        raise ValueError('Too many arguments.')
                    frequencies = None
                    array = freq_psd
                else:   # freqs and psd passed separately
                    if freqs is not None:
                        raise ValueError('Too many arguments.')
                    frequencies = freq_psd
                    array = psd
            else:
                raise ValueError('arguments not valid')
        else:   #ignore freq_psd variable
            frequencies = freqs
            array = psd

        self.psd = None
        self.freqs = None
        self.fpsd = None
        self.flogpsd = None
        self.PSD_FILTER_W = None
        self.PSD_FILTER_W_THZ = None
        self.PSD_FILTER_WF = None

        # PSD
        if array is None:
            return
        self.psd = np.array(array, dtype=float)
        self.logpsd = np.log(self.psd)
        self.logpsd_min = np.min(self.psd)

        # frequencies
        self.NFREQS = self.psd.size
        if frequencies is None:   # recompute frequencies
            self.freqs = np.linspace(0., 0.5, self.NFREQS)
        else:
            self.freqs = np.array(frequencies, dtype=float)
            if (self.freqs.size != self.NFREQS):
                raise ValueError('Number of frequencies different from PSD array size.')

        # freqs conversions to THz
        if DT_FS is not None:
            self.DT_FS = DT_FS
        self.freqs_THz = freq_red_to_THz(self.freqs, self.DT_FS)
        self.Nyquist_f_THz = self.freqs_THz[-1]
        self.DF = 0.5 / (self.NFREQS - 1)
        self.DF_THZ = freq_red_to_THz(self.DF, self.DT_FS)

    #############################################
    ###################################
    ###  COMPUTE METHODS
    ###################################
    #############################################

    def timeseries(self):
        """Return a time series (fs units)."""
        return np.arange(self.N) * self.DT_FS

    def compute_trajectory(self):
        """Compute trajectory from spectrum by IFFT."""
        if self.spectr is None:
            raise ValueError('Spectrum not defined.')
        full_spectr = np.append(self.spectr, self.spectr[-2:0:-1].conj())
        self.traj = np.real(np.fft.ifft(full_spectr))   #*np.sqrt(self.NFREQS-1)
        self.N = self.traj.size

    def compute_spectrum(self):
        """Compute spectrum from trajectory by FFT."""
        if self.traj is None:
            raise ValueError('Trajectory not defined.')
        full_spectr = np.fft.fft(self.traj)
        self.spectr = full_spectr[:self.N / 2 + 1]
        self.NFREQS = self.spectr.size
        self.DF = 0.5 / (self.NFREQS - 1)

    def compute_psd(self, PSD_FILTER_W=None, freq_units='THz', method='trajectory', DT_FS=None, normalize=False):
        # overridden in HeatCurrent (will call, at the end, this method)
        """
        Compute the periodogram from the trajectory or the spectrum.
        If a PSD_FILTER_W (expressed in freq_units) is known or given, the psd is also filtered.
        The PSD is multiplied by DT_FS at the end.
        """
        if DT_FS is not None:
            self.DT_FS = DT_FS
        if (method == 'trajectory'):
            if self.traj is None:
                raise ValueError('Trajectory not defined.')
            self.freqs, self.psdALL = periodogram(self.traj, detrend=None, axis=0)
            self.psd = np.mean(self.psdALL, axis=1)
            self.psd[1:-1] = self.psd[1:-1] * 0.5
            self.psd *= self.DT_FS
            self.NFREQS = self.freqs.size
            self.DF = 0.5 / (self.NFREQS - 1)
            self.DF_THZ = freq_red_to_THz(self.DF, self.DT_FS)
        elif (method == 'spectrum'):
            if self.spectr is None:
                raise ValueError('Spectrum not defined.')
            self.psd = self.DT_FS * np.abs(self.spectr)**2 / (2 * (self.NFREQS - 1))
            self.freqs = np.linspace(0., 0.5, self.NFREQS)
        else:
            raise KeyError('method not understood')

        self.freqs_THz = self.freqs / self.DT_FS * 1000.
        self.Nyquist_f_THz = self.freqs_THz[-1]
        if normalize:
            self.psd = self.psd / np.trapz(self.psd) / self.N / self.DT_FS
        self.logpsd = np.log(self.psd)
        self.psd_min = np.min(self.psd)
        self.psd_power = np.trapz(self.psd)   # one-side PSD power

        # (re)compute filtered psd, if a window has been defined
        if (PSD_FILTER_W is not None) or (self.PSD_FILTER_W is not None):
            self.filter_psd(PSD_FILTER_W, freq_units)

    def filter_psd(self, PSD_FILTER_W=None, freq_units='THz', window_type='rectangular', logpsd_filter_type=1):
        """
        Filter the periodogram with the given PSD_FILTER_W [freq_units].
          - PSD_FILTER_W  PSD filter window [freq_units]
          - freq_units    frequency units   ['THz', 'red' (default)]
          - window_type   filtering window type ['rectangular']
        """
        if self.psd is None:
            raise ValueError('Periodogram is not defined.')
        if PSD_FILTER_W is not None:
            if freq_units in ('THz', 'thz'):
                self.PSD_FILTER_W_THZ = PSD_FILTER_W
                self.PSD_FILTER_W = freq_THz_to_red(PSD_FILTER_W, self.DT_FS)
            elif (freq_units == 'red'):
                self.PSD_FILTER_W = PSD_FILTER_W
                self.PSD_FILTER_W_THZ = freq_red_to_THz(PSD_FILTER_W, self.DT_FS)
            else:
                raise ValueError('Freq units not valid.')
        else:
            pass   # try to use the internal value
        if self.PSD_FILTER_W is not None:
            self.PSD_FILTER_WF = int(round(self.PSD_FILTER_W * self.NFREQS * 2.))
        else:
            raise ValueError('Filter window width not defined.')

        if (window_type == 'rectangular'):
            self.fpsd = runavefilter(self.psd, self.PSD_FILTER_WF)

            # try to filter the other currents (if present)
            if self.cospectrum is not None:
                self.fcospectrum = []
                for i in range(self.cospectrum.shape[0]):
                    self.fcospectrum.append([])
                    for j in range(self.cospectrum.shape[1]):
                        ffpsd = runavefilter(self.cospectrum[i, j], self.PSD_FILTER_WF)
                        self.fcospectrum[i].append(ffpsd / self.L)
                self.fcospectrum = np.asarray(self.fcospectrum)

            # filter log-psd
            if (logpsd_filter_type == 1):
                self.flogpsd = runavefilter(self.logpsd, self.PSD_FILTER_WF)
            else:
                self.flogpsd = np.log(self.fpsd)
        else:
            raise KeyError('Window type unknown.')

    def compute_acf(self, NLAGS=None):
        """Computes the autocovariance function of the trajectory."""
        if NLAGS is not None:
            self.NLAGS = NLAGS
        else:
            self.NLAGS = self.N
        self.acf = np.zeros((self.NLAGS, self.N_COMPONENTS))
        for d in range(self.N_COMPONENTS):
            self.acf[:, d] = acovf(self.traj[:, d], unbiased=True, fft=True)[:NLAGS]
        self.acfm = np.mean(self.acf, axis=1)   # average acf

    def compute_gkintegral(self):
        """Compute the integral of the autocovariance function."""
        if self.acf is None:
            raise RuntimeError('Autocovariance is not defined.')
        self.tau = integrate_acf(self.acf)
        self.taum = np.mean(self.tau, axis=1)   # average tau

    # TODO: move this to (Heat)Current
    def compute_kappa_multi(self, others, PSD_FILTER_W=None, freq_units='red', method='trajectory', DT_FS=None,
                            normalize=False, call_other=True):   # yapf: disable
        # called by HeatCurrent
        """
        For multi-component (many current) systems: compute the cospectrum matrix and the transport coefficient.
        The results have almost the same statistical properties. The chi-square distribution has ndf = n - l + 1,
        where n is the number of time series, l is the number of currents (stored in self.ndf_chi).
         ! NOTICE: if n < l this will not work.

        In this routine the mean over the number of temporal series is already multiplied by the correct factor (the
        transport coefficient will be obtained by multiplying the result by 0.5, as in the one-component case).
        The output arrays are the same as in the one-component case.
        If a PSD_FILTER_W is known or given, the psd is also filtered.
        The elements of the matrix are multiplied by DT_FS at the end.
        others is a list of other currents, i.e. MDSample objects. For example, in the case of 4 currents, of which
        j is the energy current and j1, j2, j3 are mass currents (MDSample objects):

           j.compute_kappa_multi([j1,j2,j3], PSD_FILTER_W, freq_units)
        """
        # check if others is an array
        if not isinstance(others, (list, tuple, np.ndarray)):
            others = [others]
        N_CURRENTS = len(others)
        if DT_FS is not None:
            self.DT_FS = DT_FS

        if (method == 'trajectory'):
            if self.traj is None:
                raise ValueError('Trajectory not defined.')
            self.spectrALL = np.fft.rfft(self.traj, axis=0)
            self.NFREQS = self.spectrALL.shape[0]
            self.freqs = np.linspace(0., 0.5, self.NFREQS)
            self.DF = 0.5 / (self.NFREQS - 1)
            self.DF_THZ = freq_red_to_THz(self.DF, self.DT_FS)
        else:
            raise KeyError('method not understood')
        self.freqs_THz = self.freqs / self.DT_FS * 1000.
        self.Nyquist_f_THz = self.freqs_THz[-1]

        # calculate the same thing on the other trajectory
        if (call_other):
            for other in others:   # call other.compute_kappa_multi (MDsample method)
                MDSample.compute_kappa_multi(other, [self], PSD_FILTER_W, freq_units, method, self.DT_FS, normalize,
                                             False)
        else:
            return

        # define the cospectrum matrix. Its shape is (2, 2, NFREQS,n_spatial_dim)
        #  [  self.spectrALL*self.spectrALL.conj()     self.spectrALL*other.spectrALL.conj() ]
        #  [ other.spectrALL*self.spectrALL.conj()    other.spectrALL*other.spectrALL.conj() ]
        other_spectrALL = []
        for other in others:
            other_spectrALL.append(other.spectrALL)

        # compute the matrix defined by the outer product of only the first indexes of the two arrays
        covarALL = self.DT_FS / (2.*(self.NFREQS - 1.)) *\
                    np.einsum('a...,b...->ab...', np.array([self.spectrALL] + other_spectrALL),
                                                  np.array([self.spectrALL] + other_spectrALL).conj())

        # number of degrees of freedom of the chi-square distribution of the psd
        ndf_chi = covarALL.shape[3] - len(other_spectrALL)

        # compute the sum over the last axis (equivalent cartesian components):
        self.L = covarALL.shape[3]   ## isn't is == to N_CURRENTS?
        self.cospectrum = covarALL.sum(axis=3)

        # compute the element 1/"(0,0) of the inverse" (aka the transport coefficient)
        # the diagonal elements of the inverse have very convenient statistical properties
        multi_psd = (np.linalg.inv(self.cospectrum.transpose((2, 0, 1)))[:, 0, 0]**-1).real / ndf_chi

        if normalize:
            multi_psd = multi_psd / np.trapz(multi_psd) / self.N / self.DT_FS

        self.ndf_chi = ndf_chi
        self.psd = multi_psd
        self.logpsd = np.log(self.psd)
        self.psd_min = np.min(self.psd)
        self.psd_power = np.trapz(self.psd)   # one-side PSD power
        if (PSD_FILTER_W is not None) or (self.PSD_FILTER_W is not None):
            self.filter_psd(PSD_FILTER_W, freq_units)

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
        plot         = plot the PSD [False]
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
        return resample_timeseries(self, TSKIP, fstar_THz, FILTER_W, plot, PSD_FILTER_W, freq_units, FIGSIZE, verbose)


################################################################################

# set the default plotter of this class
MDSample.set_plotter()

################################################################################
