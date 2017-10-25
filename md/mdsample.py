
import matplotlib.pyplot as plt
from .md.tools import integrate_acf, runavefilter

class MDSample(object):
    """An MDSample object contains all the information that represent a
    single and unique Molecular Dynamics sample. 
    For example it may contain:
     - a trajectory (referred as any 1D time series 
       in real space)
     - its spectrum (the Fourier transform)
     - its Power Spectral Density, also referred as periodogram
     - ...
    All the information contained in this object should be always consistent,
    i.e. it should represent a single MD sample. Any operation that alters 
    any of the sample's properties should create a new MDSample object, in
    order to preserve the 1:1 correspondence among the sample's properties.

    An MDSample object can be initialized from any of its main properties
    (trajectory, spectrum).
    For now it is also possible to initialize a MDSample from its psd, although
    a psd does not contain enough information to regenarate a unique trajectory.

    ATTRIBUTES:
       - traj       the trajectory (any 1D real time series). 
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
       - freqsTHz   an array of frequencies, expressed in THz
       
    """
    
    def __init__(self, *args, **kwargs):
        # args -- tuple of anonymous arguments **NOT USED
        # kwargs -- dictionary of named arguments
        
        #for key, arg in kwargs.iterkeys():
        #    if (key == 'traj'):
        #        self.initialize_traj(arg)

        #check_keys = ['traj', 'spectrum']
        #for key in check_keys:
        #   arg = kwargs.get(key, None)
        #   if arg is not None:
        #     ## call correct initialize function
              
        self.traj   = kwargs.get('traj', None)
        self.spectr = kwargs.get('spectrum', None)
        #self.psd    = kwargs.get('psd', None)
        #self.fpsd   = kwargs.get('fpsd', None)
        #self.freqs  = kwargs.get('freqs', None)
        if self.traj is not None:
            self.N  = self.traj.shape[0]
            if (len(self.traj.shape) > 1):
               self.MULTI_COMPONENT = True
               self.N_COMPONENTS = self.traj.shape[1]
            else:
               self.MULTI_COMPONENT = False
               self.N_COMPONENTS = 1
        else:
            self.N  = None
        if self.spectr is not None:
            self.Nfreqs = self.spectr.size
            self.DF = 0.5 / (self.Nfreqs-1)
        #elif self.psd is not None:
        #    self.Nfreqs = self.psd.size
        else:
            self.Nfreqs = None
            self.DF = None
        self.psd = None
        self.fpsd = None
        self.logpsd = None
        self.flogpsd = None
        self.freqs = None
        self.freqsTHz = None
        self.acf = None
        self.NLAGS = None

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
    

    ###################################
    ###  INITIALIZE METHODS
    ###################################

    def initialize_traj(self, array):
        self.traj = np.array(array, dtype=float)
        self.N = self.traj.size
        if (N%2 == 1):
            raise NotImplemented('Trajectory has odd number of points.')
        if (len(self.traj.shape) > 1):
           self.MULTI_COMPONENT = True
           self.N_COMPONENTS = self.traj.shape[1]
        else:
           self.MULTI_COMPONENT = False
        return
    

    def initialize_spectrum(self, array):
        self.spectr = np.array(array, dtype=complex)
        self.Nfreqs = self.spectr.size
        self.DF = 0.5/(self.Nfreqs-1)
        return


    def initialize_psd(self, freq_psd, array=None, DT=1, DT_FS=1.0): ## modificare in (self, freq_psd=None, freqs=None, psd=None)
        if array is None:
            if (len(freq_psd) == 2):   # (freqs, psd) tuple was passed
                frequencies = freq_psd[0]
                array = freq_psd[1]
            elif (len(freq_psd) > 2):  # only psd was passed
                frequencies = None
                array = freq_psd
            else:
                raise ValueError('arguments not understood.')
        else:  # psd passed through 3rd argument (array)
            if (len(freq_psd) > 2):
                frequencies = freq_psd
            else:
                raise ValueError('Too many arguments.')
        self.psd = np.array(array, dtype=float) * DT
        self.logpsd = np.log(self.psd)
        self.psd_min = np.min(self.psd)
        self.Nfreqs = self.psd.size
        self.DF = 0.5/(self.Nfreqs-1)
        if frequencies is not None:
            self.freqs = np.array(frequencies, dtype=float)
            if (self.freqs.size != self.Nfreqs):
                raise ValueError('Number of frequencies different from PSD array size.')
        else:  # recompute frequencies
            self.freqs = np.linspace(0., 0.5, self.Nfreqs)
        self.freqsTHz = self.freqs/DT_FS*1000.
        return


    ###################################
    ###  COMPUTE METHODS
    ###################################

    def compute_trajectory(self):
        """Computes trajectory from spectrum."""
        if self.spectr is None:
            raise ValueError('Spectrum not defined.')
        full_spectr = np.append(self.spectr, self.spectr[-2:0:-1].conj())
        self.traj = np.real( ifft(full_spectr) )#*np.sqrt(self.Nfreqs-1)
        self.N = self.traj.size
        return


    def compute_spectrum(self):
        """Computes spectrum from trajectory."""
        if self.traj is None:
            raise ValueError('Trajectory not defined.')
        full_spectr = fft(self.traj)
        self.spectr = full_spectr[:self.N/2+1]
        self.Nfreqs = self.spectr.size
        self.DF = 0.5 / (self.Nfreqs-1)
        return

    
    def compute_psd(self, FILTER_WINDOW_WIDTH=None, method='trajectory', DT=1, DT_FS=1.0, average_components=True, normalize=False):
        """Computes the periodogram from the trajectory or the spectrum. 
        If a FILTER_WINDOW_WIDTH is known or given, the psd is also filtered.
        The PSD is multiplied by DT at the end."""
        if (method == 'trajectory'):
            if self.traj is None:
                raise ValueError('Trajectory not defined.')
            if self.MULTI_COMPONENT:
                self.freqs, self.psdALL = periodogram(self.traj, detrend=None, axis=0)
                self.psd = np.mean(self.psdALL, axis=1)
            else:
                self.freqs, self.psd = periodogram(self.traj, detrend=None)
            self.psd[1:-1] = self.psd[1:-1] * 0.5
            self.psd = DT * self.psd
            self.freqsTHz = self.freqs/DT_FS*1000.
            self.Nfreqs = self.freqs.size
            self.DF = 0.5 / (self.Nfreqs-1)
        elif (method == 'spectrum'):
            if self.spectr is None:
                raise ValueError('Spectrum not defined.')
            self.psd = DT * np.abs(self.spectr)**2 / (2*(self.Nfreqs - 1))
            #self.psd[1:-1] = self.psd[1:-1] * 2.0   # factor 2 from one-sided psd
            self.freqs = np.linspace(0., 0.5, self.Nfreqs)
            self.freqsTHz = self.freqs/DT_FS*1000.
        else:
            raise KeyError('method not understood')
        if normalize:
            self.psd = self.psd / np.trapz(self.psd) / self.N / DT
        self.logpsd = np.log(self.psd)
        self.psd_min = np.min(self.psd)
        if (FILTER_WINDOW_WIDTH is not None) or (self.FILTER_WINDOW_WIDTH is not None):
            self.filter_psd( FILTER_WINDOW_WIDTH )
        return


    def filter_psd(self, FILTER_WINDOW_WIDTH=None, window_type='rectangular'):
        """Filters the periodogram with the given FILTER_WINDOW_WIDTH."""
        if self.psd is None:
            raise ValueError('Periodogram is not defined.')
        if FILTER_WINDOW_WIDTH is not None:   # otherwise try to use the internal value
            self.FILTER_WINDOW_WIDTH = FILTER_WINDOW_WIDTH
        if self.FILTER_WINDOW_WIDTH is not None:
            self.FILTER_WF = int(round(self.FILTER_WINDOW_WIDTH*self.Nfreqs*2.))
        else:
            raise ValueError('Filter window width not defined.')
        if (window_type == 'rectangular'):
            self.fpsd = runavefilter(self.psd, self.FILTER_WF)
            #self.flogpsd = runavefilter(self.logpsd + EULER_GAMMA, self.FILTER_WF)
            self.flogpsd = np.log(self.fpsd)
        else:
            raise KeyError('Window type unknown.')
        return


    def compute_acf(self, NLAGS=None):
        """Computes the autocovariance function of the trajectory."""
        self.NLAGS = NLAGS
        self.acf = np.zeros((self.NLAGS+1, self.N_COMPONENTS))
        for d in xrange(self.N_COMPONENTS):
            self.acf[:,d] = acf(self.traj[:,d], nlags=NLAGS, unbiased=True, \
                                       fft=True) * np.var(self.traj[:,d])
        self.acfm = np.mean(self.acf, axis=1)  # average acf
        return


    def compute_gkintegral(self):
        """Compute the integral of the autocovariance function."""
        if self.acf is None:
            raise RuntimeError('Autocovariance is not defined.')
        self.tau = integrate_acf(self.acf)
        self.taum = np.mean(self.tau, axis=1)  # average tau
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
        plt.xticks(np.linspace(0.,0.5,11))
        plt.legend()
        return
    
    
    def plot_logpsd(self, param_dict={'label': 'log(psd)'}):
        """Plot the periodogram."""
        if self.logpsd is None:
            raise ValueError('Log-Peridogram not defined.')
        plt.plot(self.freqs, self.logpsd, **param_dict)
        plt.xlabel('f [$\omega$*DT/2$\pi$]')
        plt.xticks(np.linspace(0.,0.5,11))
        plt.legend()
        return
    
    
    def plot_fpsd(self, FILTER_WINDOW_WIDTH=None, param_dict={'label': 'f-psd'}):
        """Plot the filtered periodogram. 
        If FILTER_WINDOW_WIDTH is defined/passed a filtered psd is computed,
        otherwise the internal copy is used.
        """
        if (FILTER_WINDOW_WIDTH is not None) or (self.FILTER_WINDOW_WIDTH is not None):
            self.filter_psd( FILTER_WINDOW_WIDTH )
        if self.fpsd is None:
            raise ValueError('Filtered peridogram not defined.')
        plt.plot(self.freqs, self.fpsd, **param_dict)
        plt.xlabel('f [$\omega$*DT/2$\pi$]')
        plt.xticks(np.linspace(0.,0.5,11))
        plt.legend()
        return

        
    def plot_flogpsd(self, FILTER_WINDOW_WIDTH=None, param_dict={'label': 'f-log(psd)'}):
        """Plot the filtered periodogram. 
        If FILTER_WINDOW_WIDTH is defined/passed a filtered psd is computed,
        otherwise the internal copy is used.
        """
        if (FILTER_WINDOW_WIDTH is not None) or (self.FILTER_WINDOW_WIDTH is not None):
            self.filter_psd( FILTER_WINDOW_WIDTH )
        if self.flogpsd is None:
            raise ValueError('Filtered log-peridogram not defined.')
        plt.plot(self.freqs, self.flogpsd, **param_dict)
        plt.xlabel('f [$\omega$*DT/2$\pi$]')
        plt.xticks(np.linspace(0.,0.5,11))
        plt.legend()
        return
