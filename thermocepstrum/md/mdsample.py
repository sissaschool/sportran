import numpy as np

#import matplotlib.pyplot as plt
from thermocepstrum.utils.loadAfterPlt import plt

from .tools import integrate_acf, runavefilter
from scipy.signal import periodogram
from .acf import acovf


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
                    Nfreqs = N/2+1 normalized frequencies contained
                    in the interval [0, 1/(2N*DT)]
       - psd        the Power Spectral Density (periodogram), defined as
                              DT    N-1
                     I(f) =  ---- * SUM | x[n] * exp(-2.0J*pi*f/N)
                               N    n=0
                    with f = [0, 1/(2N*DT)]
       - N          size of traj
       - Nfreqs     number of trajectories, should be N/2+1
       - freqs      an array of frequencies, should be [0, 1/(2N*DT)]
       - freqs_THz  an array of frequencies, expressed in THz

    MEMBERS:
        self.DT_FS                  timestep in femtoseconds
        self.fpsd                   filtered periodogram
        self.flogpsd                filtered log-periodogram
        self.acf                    autocorrelation function
        self.N_COMPONENTS           number of equivalent (cartesian) components (an average over them will be computed)
        self.MULTI_COMPONENT        True if N_COMPONENTS > 1
        self.FILTER_WINDOW_WIDTH    width of the moving average filter (reduced frequency units)
        self.FILTER_WF              width of the moving average filter (number of frequencies)

    """

    def __init__(self, traj=None, spectr=None, psd=None, freqs=None, DT_FS=1.0):
        self.DT_FS = DT_FS
        self.initialize_traj(traj)
        self.initialize_spectrum(spectr)
        self.initialize_psd(freqs=freqs, psd=psd, DT_FS=DT_FS)

        self.fpsd = None
        self.flogpsd = None
        self.acf = None
        self.NLAGS = None
        self.cospectrum = None
        self.fcospectrum = None

        # other variables...
        self.FILTER_WINDOW_WIDTH = None
        self.FILTER_WF = None
        return

    def __repr__(self):
        msg = 'MDSample:\n' + \
              '  traj:   {}  steps  *  {} components\n'.format(self.N, self.N_COMPONENTS) + \
              '  spectr: {}  frequencies\n'.format(self.Nfreqs)
        if self.psd is not None:
            msg = msg + '  psd:    {}  frequencies\n'.format(self.psd.size) + \
                        '    DF =   {}  [omega*DT/2/pi]\n'.format(self.DF)
        if self.fpsd is not None:
            msg = msg + '  fpsd:   {}  frequencies\n'.format(self.fpsd.size) +\
                  '    FILTER_WINDOW_WIDTH = {} [omega*DT/2/pi]\n'.format(self.FILTER_WINDOW_WIDTH) +\
                  '    FILTER_WF           = {} frequencies\n'.format(self.FILTER_WF)
        if self.acf is not None:
            msg = msg + '  acf:    {}  lags\n'.format(self.NLAGS)
        return msg

    #############################################
    ###################################
    ###  INITIALIZE METHODS
    ###################################
    #############################################

    def initialize_traj(self, array):
        if array is not None:
            if array.shape[0] % 2 == 1:
                self.traj = np.array(array[1:], dtype=float)
                print('trajectory has an odd number of points. Removing the first one.')
            else:
                self.traj = np.array(array, dtype=float)
            self.N = self.traj.shape[0]
            if len(self.traj.shape) > 1:
                self.MULTI_COMPONENT = True
                self.N_COMPONENTS = self.traj.shape[1]
            else:
                self.MULTI_COMPONENT = False
                self.N_COMPONENTS = 1
        else:
            self.traj = None
            self.N = None
            self.N_COMPONENTS = None
        return

    def initialize_spectrum(self, array):
        if array is not None:
            self.spectr = np.array(array, dtype=complex)
            self.Nfreqs = self.spectr.size
            self.DF = 0.5 / (self.Nfreqs - 1)
        else:
            self.spectr = None
            self.Nfreqs = None
            self.DF = None
        return

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
        else:   #ignore freq)psd variable
            frequencies = freqs
            array = psd

        # PSD
        if array is None:
            self.psd = None
            self.freqs = None
            return
        self.psd = np.array(array, dtype=float)
        self.logpsd = np.log(self.psd)
        self.logpsd_min = np.min(self.psd)

        # frequencies
        self.Nfreqs = self.psd.size
        if frequencies is None:   # recompute frequencies
            self.freqs = np.linspace(0., 0.5, self.Nfreqs)
        else:
            self.freqs = np.array(frequencies, dtype=float)
            if (self.freqs.size != self.Nfreqs):
                raise ValueError('Number of frequencies different from PSD array size.')

        # freqs conversions to THz
        if DT_FS is not None:
            self.DT_FS = DT_FS
        self.freqs_THz = self.freqs / self.DT_FS * 1000.
        self.Nyquist_f_THz = self.freqs_THz[-1]
        self.DF = 0.5 / (self.Nfreqs - 1)
        self.DF_THz = self.DF / self.DT_FS * 1000.
        return

    #############################################
    ###################################
    ###  COMPUTE METHODS
    ###################################
    #############################################

    def timeseries(self):
        return np.arange(self.N) * self.DT_FS

    def compute_trajectory(self):
        """Computes trajectory from spectrum."""
        if self.spectr is None:
            raise ValueError('Spectrum not defined.')
        full_spectr = np.append(self.spectr, self.spectr[-2:0:-1].conj())
        self.traj = np.real(np.fft.ifft(full_spectr))   #*np.sqrt(self.Nfreqs-1)
        self.N = self.traj.size
        return

    def compute_spectrum(self):
        """Computes spectrum from trajectory."""
        if self.traj is None:
            raise ValueError('Trajectory not defined.')
        full_spectr = np.fft.fft(self.traj)
        self.spectr = full_spectr[:self.N / 2 + 1]
        self.Nfreqs = self.spectr.size
        self.DF = 0.5 / (self.Nfreqs - 1)
        return

    #overridden in HeatCurrent (will call, at the end, this method)
    def compute_psd(self, FILTER_WINDOW_WIDTH=None, method='trajectory', DT_FS=None, average_components=True,
                    normalize=False):   # yapf: disable
        """
        Compute the periodogram from the trajectory or the spectrum.
        If a FILTER_WINDOW_WIDTH (reduced frequency units) is known or given, the psd is also filtered.
        The PSD is multiplied by DT_FS at the end.
        """

        if DT_FS is not None:
            self.DT_FS = DT_FS
        if (method == 'trajectory'):
            if self.traj is None:
                raise ValueError('Trajectory not defined.')
            if self.MULTI_COMPONENT:
                self.freqs, self.psdALL = periodogram(self.traj, detrend=None, axis=0)
                self.psd = np.mean(self.psdALL, axis=1)
            else:
                self.freqs, self.psd = periodogram(self.traj, detrend=None)
            self.psd[1:-1] = self.psd[1:-1] * 0.5
            self.psd = self.DT_FS * self.psd
            self.Nfreqs = self.freqs.size
            self.DF = 0.5 / (self.Nfreqs - 1)
        elif (method == 'spectrum'):
            if self.spectr is None:
                raise ValueError('Spectrum not defined.')
            self.psd = self.DT_FS * np.abs(self.spectr)**2 / (2 * (self.Nfreqs - 1))
            #self.psd[1:-1] = self.psd[1:-1] * 2.0   # factor 2 from one-sided psd
            self.freqs = np.linspace(0., 0.5, self.Nfreqs)
        else:
            raise KeyError('method not understood')
        self.freqs_THz = self.freqs / self.DT_FS * 1000.
        self.Nyquist_f_THz = self.freqs_THz[-1]
        if normalize:
            self.psd = self.psd / np.trapz(self.psd) / self.N / self.DT_FS
        self.logpsd = np.log(self.psd)
        self.psd_min = np.min(self.psd)
        self.psd_power = np.trapz(self.psd)   # one-side PSD power
        if (FILTER_WINDOW_WIDTH is not None) or (self.FILTER_WINDOW_WIDTH is not None):
            self.filter_psd(FILTER_WINDOW_WIDTH)
        return

    def filter_psd(self, FILTER_WINDOW_WIDTH=None, window_type='rectangular', logpsd_filter_type=1):
        """Filters the periodogram with the given FILTER_WINDOW_WIDTH [freq units]."""
        if self.psd is None:
            raise ValueError('Periodogram is not defined.')
        if FILTER_WINDOW_WIDTH is not None:   # otherwise try to use the internal value
            self.FILTER_WINDOW_WIDTH = FILTER_WINDOW_WIDTH
        if self.FILTER_WINDOW_WIDTH is not None:
            self.FILTER_WF = int(round(self.FILTER_WINDOW_WIDTH * self.Nfreqs * 2.))
        else:
            raise ValueError('Filter window width not defined.')
        if (window_type == 'rectangular'):
            self.fpsd = runavefilter(self.psd, self.FILTER_WF)
            if self.cospectrum is not None:   # try to filter the other currents (if present)
                self.fcospectrum = []
                for i in range(self.cospectrum.shape[0]):
                    self.fcospectrum.append([])
                    for j in range(self.cospectrum.shape[1]):
                        ffpsd = runavefilter(self.cospectrum[i, j], self.FILTER_WF)
                        self.fcospectrum[i].append(ffpsd / self.L)
                self.fcospectrum = np.asarray(self.fcospectrum)
            if logpsd_filter_type == 1:
                self.flogpsd = runavefilter(self.logpsd, self.FILTER_WF)
            else:
                self.flogpsd = np.log(self.fpsd)
        else:
            raise KeyError('Window type unknown.')
        return

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
        return

    def compute_gkintegral(self):
        """Compute the integral of the autocovariance function."""
        if self.acf is None:
            raise RuntimeError('Autocovariance is not defined.')
        self.tau = integrate_acf(self.acf)
        self.taum = np.mean(self.tau, axis=1)   # average tau
        return

    # this is called by HeatCurrent.
    def compute_kappa_multi(self, others, FILTER_WINDOW_WIDTH=None, method='trajectory', DT_FS=None,
                            average_components=True, normalize=False, call_other=True):   # yapf: disable
        """
        For multi-component (many current) systems: compute the cospectrum matrix and the transport coefficient.
        The results have almost the same statistical properties. The chi-square distribution has ndf = n - l + 1,
        where n is the number of time series, l is the number of currents (stored in self.ndf_chi).
         ! NOTICE: if n < l this will not work.

        In this routine the mean over the number of temporal series is already multiplied by the correct factor (the
        transport coefficient will be obtained by multiplying the result by 0.5, as in the one-component case).
        The output arrays are the same as in the one-component case.
        If a FILTER_WINDOW_WIDTH is known or given, the psd is also filtered.
        The elements of the matrix are multiplied by DT_FS at the end.
        others is a list of other currents, i.e. MDSample objects. For example, in the case of 4 currents, of which
        j is the energy current and j1, j2, j3 are mass currents (MDSample objects):

           j.compute_kappa_multi(others=[j1,j2,j3], FILTER_WINDOW_WIDTH=FILTER_WINDOW_WIDTH)
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
            self.Nfreqs = self.spectrALL.shape[0]
            self.freqs = np.linspace(0., 0.5, self.Nfreqs)
            self.DF = 0.5 / (self.Nfreqs - 1)
        else:
            raise KeyError('method not understood')
        self.freqs_THz = self.freqs / self.DT_FS * 1000.
        self.Nyquist_f_THz = self.freqs_THz[-1]

        # calculate the same thing on the other trajectory
        if (call_other):
            for other in others:   # call other.compute_kappa_multi (MDsample method)
                MDSample.compute_kappa_multi(other, [self], FILTER_WINDOW_WIDTH, method, self.DT_FS, average_components,
                                             normalize, False)
        else:
            return

        # define the cospectrum matrix. Its shape is (2, 2, Nfreqs,n_spatial_dim)
        #  [  self.spectrALL*self.spectrALL.conj()     self.spectrALL*other.spectrALL.conj() ]
        #  [ other.spectrALL*self.spectrALL.conj()    other.spectrALL*other.spectrALL.conj() ]
        other_spectrALL = []
        for other in others:
            other_spectrALL.append(other.spectrALL)

        # compute the matrix defined by the outer product of only the first indexes of the two arrays
        covarALL = self.DT_FS / (2.*(self.Nfreqs - 1.)) *\
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
        if (FILTER_WINDOW_WIDTH is not None) or (self.FILTER_WINDOW_WIDTH is not None):
            self.filter_psd(FILTER_WINDOW_WIDTH)
        return

    ###################################
    ###  PLOT METHODS
    ###################################
    # customized line properties can be passed through the param_dict dictionary
    # see kwargs at: http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot

    def plot_traj(self, param_dict={'label': 'traj'}):
        """Plot the time series."""
        if self.traj is None:
            raise ValueError('Trajectory not defined.')
        plt.plot(self.traj, **param_dict)
        plt.xlabel('t')
        plt.grid()
        plt.legend()
        return

    def plot_psd(self, param_dict={'label': 'psd'}):
        """Plot the periodogram."""
        if self.psd is None:
            raise ValueError('Peridogram not defined.')
        plt.plot(self.freqs, self.psd, **param_dict)
        plt.xlabel('f [$\omega$*DT/2$\pi$]')
        plt.xticks(np.linspace(0., 0.5, 11))
        plt.legend()
        return

    def plot_logpsd(self, param_dict={'label': 'log(psd)'}):
        """Plot the periodogram."""
        if self.logpsd is None:
            raise ValueError('Log-Peridogram not defined.')
        plt.plot(self.freqs, self.logpsd, **param_dict)
        plt.xlabel('f [$\omega$*DT/2$\pi$]')
        plt.xticks(np.linspace(0., 0.5, 11))
        plt.legend()
        return

    def plot_fpsd(self, FILTER_WINDOW_WIDTH=None, param_dict={'label': 'f-psd'}):
        """Plot the filtered periodogram.
        If FILTER_WINDOW_WIDTH is defined/passed a filtered psd is computed,
        otherwise the internal copy is used.
        """
        if (FILTER_WINDOW_WIDTH is not None) or (self.FILTER_WINDOW_WIDTH is not None):
            self.filter_psd(FILTER_WINDOW_WIDTH)
        if self.fpsd is None:
            raise ValueError('Filtered peridogram not defined.')
        plt.plot(self.freqs, self.fpsd, **param_dict)
        plt.xlabel('f [$\omega$*DT/2$\pi$]')
        plt.xticks(np.linspace(0., 0.5, 11))
        plt.legend()
        return

    def plot_flogpsd(self, FILTER_WINDOW_WIDTH=None, param_dict={'label': 'f-log(psd)'}):
        """Plot the filtered periodogram.
        If FILTER_WINDOW_WIDTH is defined/passed a filtered psd is computed,
        otherwise the internal copy is used.
        """
        if (FILTER_WINDOW_WIDTH is not None) or (self.FILTER_WINDOW_WIDTH is not None):
            self.filter_psd(FILTER_WINDOW_WIDTH)
        if self.flogpsd is None:
            raise ValueError('Filtered log-peridogram not defined.')
        plt.plot(self.freqs, self.flogpsd, **param_dict)
        plt.xlabel('f [$\omega$*DT/2$\pi$]')
        plt.xticks(np.linspace(0., 0.5, 11))
        plt.legend()
        return
