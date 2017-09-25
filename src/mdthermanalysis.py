################################################################################
###
###   MDThermAnalysis - v0.3.10  - Sep 22, 2017
###
################################################################################
###
###  a package to analyze thermal currents from Equilibrium Molecular Dynamics
###  simulations
### 
################################################################################

import numpy as np
from scipy.fftpack import fft, ifft, dct
from scipy.signal import periodogram, lfilter
import matplotlib.pyplot as plt
from statsmodels.tsa.api import acf

EULER_GAMMA = 0.57721566490153286060651209008240243104215933593992  # Euler-Mascheroni constant

################################################################################

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
        return


    def compute_gkintegral(self):
        """Compute the integral of the autocovariance function."""
        if self.acf is None:
            raise RuntimeError('Autocovariance is not defined.')
        self.tau = integrate_acf(self.acf)
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

################################################################################
################################################################################

def integrate_acf(acf):
    """Returns the integral function of f, i.e. its integral as a function of 
    the upper integration limit. Supports multi-component arrays."""
    N = acf.shape[0]
    tau = np.zeros(acf.shape)
    for i in range(1, N):
        tau[i] = tau[i-1] + 0.5*acf[i-1] + 0.5*acf[i]
    return tau


def filter_and_sample( y_big, W, DT, window='rectangular', even_NSTEPS=True, detrend=False, drop_first=True ):
    """Filter signal with moving average window of width W and sample it with time step DT."""
    
    if ( W > 1 ):
        if (window == 'rectangular'):
            y_f = lfilter( (1./W)*np.ones(W), 1., y_big, axis=0 );
        else:
            raise NotImplemented('Not implemented window type.')
        
        # drop first W steps (initial conditions)
        if drop_first:
            y = y_f[ (W-1)::DT ];
        else:
            y = y_f[ ::DT ];
    else:
        y = y_big[ ::DT ];
    
    # remove the mean
    if detrend:
        y = y - np.mean(y, axis=0)
    
    # keeps an even number of points
    if even_NSTEPS:
        if (y.shape[0]%2 == 1):
            return y[:-1]
    return y


def runavefilter(X, WF):
    """ Computes the running average of a numpy array over WF consecutive
    elements (or WF+1 if WF is even):
            (X[i-WF/2]+...+X[i]+...+X[i+WF/2]) / WF
    assumes elf.tau_std_Kmin = np.sqrt(self.tau_var_Kmin)
                                                                        hat the array is "even" ( X[-i] = X[i] ) and anti-periodic
    (X[(N-1)+i]=X[(N-1)-i]), like a ONE-SIDED PSD."""
    
    if (WF%2 == 0):
        WF = WF + 1
    W = int(WF/2)
    
    Y = np.concatenate( (X[W:0:-1], X, X[-2:-W-2:-1]) )
    return np.convolve( Y, np.array([1.0/WF]*WF), 'valid' )


def generate_empirical_spectrum(psd):
    """Add noise to a periodogram and generate a complex spectrum."""
    ### Probabilmente ci vuole un fattore 1/N, in maniera da tirar via il fattore N 
    ### in compute_trajectory, cosi' diventa consistente con compute_spectrum
    spectr = np.random.normal(0., np.sqrt(0.5*psd*(psd.size-1)*2)) + \
               1.0J*np.random.normal(0., np.sqrt(0.5*psd*(psd.size-1)*2))
    spectr[0]  = spectr[0].real
    spectr[-1] = spectr[-1].real
    return spectr


################################################################################
################################################################################

class AR_Model(object):
    """An AR_Model defines an AutoRegressive process of order P.
    It is possible to fit a time series to an AR(P) model of order P.
    
    ATTRIBUTES:
       - P          the order of the AutoRegressive process.
       - phi        the P parameters of the process.
       - sigma2     the variance of the innovations.

    """

    def __init__(self, *args, **kwargs):
        # args -- tuple of anonymous arguments **NOT USED
        # kwargs -- dictionary of named arguments
        
        if (len(args) == 1):
            self.P = args[0]
        else:
            self.P   = kwargs.get('order', None)
        self.phi = kwargs.get('phi', None)
        self.sigma2 = kwargs.get('sigma2', None)
        if self.phi is not None:
            if self.P is None:
                self.P = self.phi.size
            elif self.phi.size != self.P:
                raise ValueError('phi array size is different from the given P.')
        self.phi_std = None
        self.sigma2_std = None
        self.cov = None
        return
    
    
    def __repr__(self):
        msg = 'AR({}) model:\n'.format(self.P) + \
              '  phi = {}\n'.format(self.phi) + \
              '  sigma2 = {}\n'.format(self.sigma2)
        return msg
    
    
    def fit(self, traj, order=None):
        if order is None:
            order = self.P
            if order is None:
                raise ValueError('P is not defined')
        #if self.traj is None:
        #    raise ValueError('Trajectory not loaded')
        self.phi, self.sigma2, self.cov = CSS_Solve(traj, order)
        self.phi_std = np.sqrt( np.diag(self.cov)[:-1] )
        self.sigma2_std = np.sqrt( self.cov[-1,-1] )
        return
    
    
    def compute_psd(self, NFREQS, DT=1):
        """Computes the theoretical periodogram of the AR(p) model."""
        if self.phi is None:
            raise ValueError('AR phi coefficients not defined.')
        if self.sigma2 is None:
            raise ValueError('AR sigma2 not defined.')
        return ar_psd(self.phi, self.sigma2, NFREQS, DT)

    
    def compute_tau(self, DT=1):
        """Computes the theoretical zero frequency of the periodogram."""
        if self.phi is None:
            raise ValueError('AR phi coefficients not defined.')
        if self.sigma2 is None:
            raise ValueError('AR sigma2 not defined.')
        return ar_tau(self.phi, self.sigma2, self.cov, DT)
    
    
    def generate_trajectory(self, N):
        """Generates an AR(p) trajectory of length N."""
        if self.phi is None:
            raise ValueError('AR phi coefficients not defined.')
        if self.sigma2 is None:
            raise ValueError('AR sigma2 not defined.')
        traj = np.zeros(N + self.P)
        noise_std = np.sqrt(self.sigma2)
        for t in range(self.P, N + self.P):
            x_t = np.random.normal(0., noise_std)
            for i in range(self.P):
                x_t = x_t + traj[t-i-1]*self.phi[i]
            traj[t] = x_t
        return traj[self.P:]

################################################################################
################################################################################

def CSS_Solve(y, P):
    """Solve Conditionate Sum-of-Squares.

    INPUT:    y  : time series
              P  : AR order
    RETURNS:  phi    : AR coefficients
              sigma2 : variance of innovations
              V      : asymptotic covariance matrix of parameters"""
        
    from scipy.linalg import solve, inv

    RUN_TIME = y.size

    # compute AR phi
    a = np.zeros((P,P))
    for i in range(P):
        for j in range(i,P):
            c = 0.;
            for t in range(P,RUN_TIME):
                c += y[t-(i+1)]*y[t-(j+1)]/(RUN_TIME-P)
            a[i,j] = a[j,i] = c
    
    b = np.zeros(P)
    for i in range(P):
        c = 0.;
        for t in range(P,RUN_TIME):
            c += y[t-(i+1)]*y[t]/(RUN_TIME-P)
        b[i] = c
    
    phi = solve( a, b, sym_pos=True )
    
    # compute sigma2
    sigma2 = 0.
    for t in range(P,RUN_TIME):
        sigma2 += ( y[t] - np.dot(y[t-P:t],phi[::-1]) )**2
        #res = y[t]
        #for j in range(P):
        #    res -= phi[j]*y[t-(j+1)]
        #sigma2 += res**2
    sigma2 *= 1./(RUN_TIME-P)
        
    # compute asymptotic covariance matrix of parameters
    yy = y[P-1:-1]
    for i in range(P-1):
        yy = np.row_stack( (yy, y[P-2-i:-2-i]) )
        
    V = np.zeros((P+1,P+1))
    V[:P,:P] = inv( np.cov( yy ) )*sigma2/(RUN_TIME-P-1)
    V[P,P]   = 2.*sigma2*sigma2/(RUN_TIME-P-1)
    return phi, sigma2, V


def ar_psd(AR_phi, AR_sigma2, NFREQS, DT=1):
    """Compute psd of an AR(P) process. BE CAREFUL WITH NORMALIZATION IF DT!=1"""
    P = len(AR_phi)
    freqs  = np.linspace( 0., 0.5/DT, NFREQS+1 )
    AR_psd = np.zeros( freqs.size )
    for i in range( freqs.size ):
        phiz = np.sum( AR_phi * np.exp( -2.0j * np.pi * freqs[i] * np.arange(1, P + 1) ) )
        #AR_psd[i] = 2. * DT * AR_sigma2 / np.abs( 1. - phiz )**2  # factor 2 comes from one-sided
        AR_psd[i] = DT * AR_sigma2 / np.abs( 1. - phiz )**2
    return freqs/DT, AR_psd


def ar_tau(AR_phi, AR_sigma2, AR_phi_cov=None, DT=1, RUN_TIME=0):
    """Compute tau of an AR(P) process."""
    
    P = AR_phi.size;
    
    if AR_phi_cov is not None:
        # check if covariance matrix contains cov(sigma2,_)
        if ( AR_phi_cov.shape == (P+1,P+1) ):
            AR_cov = AR_phi_cov
        elif AR_phi_cov.shape == (P,P):
            if ( RUN_TIME == 0 ):
                raise ValueError('You should pass RUN_TIME to AR_tau function.')
            # build total covariance matrix (phi+sigma2)
            AR_cov          = np.zeros( (P+1,P+1) )
            AR_cov[:-1,:-1] = AR_phi_cov;
            AR_cov[-1,-1]   = 2. * AR_sigma2**2 / (RUN_TIME-P-1);
        else:
            raise ValueError('AR_phi_cov matrix has wrong dimensions.')

    # compute model tau from S(f=0)
    phiz  = np.sum( AR_phi );
    AR_tau   = 0.5 * DT * AR_sigma2 / np.abs( 1. - phiz )**2;
    
    if AR_phi_cov is not None:
        grad_tau = DT * np.append( AR_sigma2/phiz**3*np.ones(P), 0.5/phiz**2 );
        AR_tau_stderr = np.sqrt( grad_tau.dot(AR_cov).dot(grad_tau) );
        return AR_tau, AR_tau_stderr
    else:
        return AR_tau


################################################################################
################################################################################

class CosFilter(object):
    
    def __init__(self, samplelogpsd, ck_theory_var=None, psd_theory_mean=None, aic_type='aic', Kmin_corrfactor=1.0, normalization=1.0):

        NF = samplelogpsd.size
        N = 2*(NF-1)

        if psd_theory_mean is None:
            # by default the THEORETICAL means are the one component ones:
            # ck THEORY mean:
            #    - EULER_GAMMA - log(2)   for k = {0, N/2}
            #    - EULER_GAMMA            otherwise
            self.logpsd_THEORY_mean = - EULER_GAMMA * np.ones(NF)
            self.logpsd_THEORY_mean[0]  = - EULER_GAMMA - np.log(2)
            self.logpsd_THEORY_mean[-1] = - EULER_GAMMA - np.log(2)
        else:
            self.logpsd_THEORY_mean = psd_theory_mean

        # subtract the mean of the distribution
        self.samplelogpsd = samplelogpsd - self.logpsd_THEORY_mean
        self.logpsdK = dct_coefficients(self.samplelogpsd, normalization)
        if (aic_type == 'aic'):
            self.aic = dct_AIC(self.logpsdK, ck_theory_var)
        elif (aic_type == 'aicc'):
            self.aic = dct_AICc(self.logpsdK, ck_theory_var)
        else:
            raise ValueError('AIC type not valid.')
        self.aic_min = np.min(self.aic)
        self.aic_Kmin = int(round(np.argmin(self.aic) * Kmin_corrfactor))
        if (self.aic_Kmin >= NF):
            print "! Warning:  aic_Kmin ({:}) is out of range.".format(self.aic_Kmin)
        
        if ck_theory_var is None:
            # by default the THEORETICAL variances are the one component ones:
            # ck THEORY variances:
            #    (pi^2)/3/N   for k = {0, N/2}
            #    (pi^2)/6/N   otherwise
            self.logpsdK_THEORY_var = 1./N * np.concatenate(([np.pi**2/3], [np.pi**2/6.]*(NF-2), [np.pi**2/3]))
            self.logpsdK_THEORY_std = np.sqrt(self.logpsdK_THEORY_var)
            # logtau THEORY variances:  (we assume to be summing ck up to K, included)
            #    (pi^2)/3/N*(2*K+1)   for K = {0, N/2-1}
            #    (pi^2)/3             for K = N/2
            self.logtau_THEORY_var = 1./N * np.concatenate((np.pi**2/3.*(2*np.arange(NF-1)+1), [np.pi**2/3.*N]))
            self.logtau_THEORY_std = np.sqrt(self.logtau_THEORY_var)
        else:
            self.logpsdK_THEORY_var = ck_theory_var
            self.logpsdK_THEORY_std = np.sqrt(self.logpsdK_THEORY_var)
            self.logtau_THEORY_var = np.zeros(NF)
            self.logtau_THEORY_var[0] = self.logpsdK_THEORY_var[0]
            for K in xrange(1, NF-1):
                self.logtau_THEORY_var[K] = self.logtau_THEORY_var[K-1] + 4.*self.logpsdK_THEORY_var[K]
            self.logtau_THEORY_var[-1] = self.logtau_THEORY_var[-2] + self.logpsdK_THEORY_var[-1]
            self.logtau_THEORY_std = np.sqrt(self.logtau_THEORY_var)
        return
    
    
    def scan_filter_tau(self, K_PSD=None, correct_mean=True):
        """Computes the tau as a function of the cutoff K for the CosFilter.
        Also computes psd and logpsd for the given K_PSD cutoff (or otherwise the aic_Kmin is used)."""
        if K_PSD is None:
            self.K_PSD = self.aic_Kmin
        else:
            self.K_PSD = K_PSD

        # COS-filter analysis with frequency cutoff K
        self.logtau = dct_filter_tau(self.samplelogpsd)
        self.logpsd = dct_filter_psd(self.samplelogpsd, self.K_PSD) # usually is log(psd)@aic_Kmin
        self.psd = np.exp(self.logpsd)
        self.tau = np.exp(self.logtau)

        if (self.aic_Kmin < self.samplelogpsd.size):
            self.logtau_Kmin     = self.logtau[self.aic_Kmin]
            self.logtau_var_Kmin = self.logtau_THEORY_var[self.aic_Kmin]
            self.logtau_std_Kmin = np.sqrt(self.logtau_var_Kmin)
            self.tau_Kmin     = self.tau[self.aic_Kmin]
            self.tau_var_Kmin = self.tau_Kmin * self.logtau_var_Kmin
            self.tau_std_Kmin = np.sqrt(self.tau_var_Kmin)
        else:
            self.logtau_Kmin     = np.NaN
            self.logtau_var_Kmin = np.NaN
            self.logtau_std_Kmin = np.NaN
            self.tau_Kmin     = np.NaN
            self.tau_var_Kmin = np.NaN
            self.tau_std_Kmin = np.NaN

        if correct_mean:
            self.logpsd = self.logpsd + self.logpsd_THEORY_mean
            self.logtau = self.logtau + self.logpsd_THEORY_mean[0]
            self.logtau_Kmin = self.logtau_Kmin + self.logpsd_THEORY_mean[0]
        return
    
    
    def scan_filter_psd(self, K_LIST, correct_mean=True):
        """Computes the psd as a function of the cutoff K for the CosFilter.
        Repeats the procedure for all the cutoffs in the K_LIST."""
        self.K_LIST = K_LIST
        self.logpsd_K_LIST  = np.zeros((self.samplelogpsd.size, len(self.K_LIST)))
        self.psd_K_LIST     = np.zeros((self.samplelogpsd.size, len(self.K_LIST)))
        self.logtau_K_LIST  = np.zeros(len(self.K_LIST))   # DEFINED AS log(PSD[0]), no factor 0.5 or 0.25
        self.tau_K_LIST     = np.zeros(len(self.K_LIST))
        
        for k, K in enumerate(self.K_LIST):
            # COS-filter analysis with frequency cutoff K
            self.logpsd_K_LIST[:,k] = dct_filter_psd(self.samplelogpsd, K)
            self.logtau_K_LIST[k]   = self.logpsd_K_LIST[0,k]

            self.psd_K_LIST[:,k]  = np.exp(self.logpsd_K_LIST[:,k])
            self.tau_K_LIST[k]    = np.exp(self.logtau_K_LIST[k])

            if correct_mean:
                self.logpsd_K_LIST[:,k] = self.logpsd_K_LIST[:,k] + self.logpsd_THEORY_mean
                self.logtau_K_LIST[k]   = self.logtau_K_LIST[k] + self.logpsd_THEORY_mean[0]
        return

    
    def compute_p_aic(self, method='ba'):
        """Define a weight distribution from the AIC, according to a method."""
        NF = self.samplelogpsd.size
        self.p_aic = produce_p(self.aic, method)
        self.p_aic_Kave, self.p_aic_Kstd = grid_statistics(np.arange(NF), self.p_aic)
        return
    
    
    def compute_logtau_density(self, method='ba', only_stats=False, density_grid=None, grid_size=1000, correct_mean=True):
        if self.p_aic is None:
            raise ValueError('No P_AIC defined.')

        # compute statistics
        self.p_logtau_density_xave, self.p_logtau_density_xstd = \
                        grid_statistics(self.logtau, self.p_aic, self.logtau_THEORY_var + self.logtau**2)
        self.p_logtau_density_xstd2 = np.dot(self.p_aic, np.sqrt(self.logtau_THEORY_var + (self.logtau-self.p_logtau_density_xave)**2))
        ##self.p_logtau_density_xave, self.p_logtau_density_xstd = \
        ##                ta.grid_statistics(self.p_logtau_grid, self.p_logtau_density)

        # compute distribution
        if not only_stats:
            if density_grid is None:
                self.p_logtau_density, self.p_logtau_grid = produce_p_density(self.p_aic, \
                                        self.logtau_THEORY_std, self.logtau, grid_size=grid_size)
            else:
                self.p_logtau_grid = density_grid
                self.p_logtau_density = produce_p_density(self.p_aic, self.logtau_THEORY_std, \
                                            self.logtau, grid=self.p_logtau_grid)

        # tau distribution
        self.p_tau_density_xave, self.p_tau_density_xstd = logtau_to_tau(self.p_logtau_density_xave, self.logpsd_THEORY_mean[0], self.p_logtau_density_xstd)
        return

#    def optimize_cos_filter(self, thr=0.05, K_LIST=None, logtauref=None):
#        if K_LIST is not None:
#            self.K_LIST = K_LIST
#        self.scan_cos_filter_K()
#        ## find minimum cutoff K that satisfies  |log(tau) - tauref| < thr
#        if logtauref is not None:
#            self.logtauref = logtauref
#        else:
#            self.logtauref = self.logtau[-1]  # if tauref is not given, use logtau with max cutoff
#        self.optimalK_idx = len(self.K_LIST) - np.argmin(np.abs(self.logtau - self.logtauref)[::-1] <= thr)
#        if (self.optimalK_idx < len(self.K_LIST)):
#            self.optimalK = self.K_LIST[self.optimalK_idx]
#        else:
#            self.optimalK_idx = np.NaN
#            self.optimalK = np.NaN
#            print 'Warning: optimal cutoff K NOT FOUND.'
#        return

################################################################################


def logtau_to_tau(logtau, logtau_mean, logtau_var, correct_mean=True):
    if correct_mean:
        logtau = logtau - logtau_mean
        tau = np.exp(logtau)
        tau_std = np.sqrt(tau * logtau_var)
    else:
        raise NotImplemented()
    return tau, tau_std

def dct_coefficients(y, normalization=1.0):
    yk = dct(y, type=1)/(y.size-1)*0.5 #normalization
    return yk

def dct_AIC(yk, theory_var=None):
    """AIC[K] = sum_{k>K} c_k^2/theory_var + 2*(K+1)
    Assumiamo di tenere tutti i k <= K."""
    aic = np.zeros(yk.size)
    if theory_var is None:
        N = 2 * (yk.size - 1)
        theory_var1 = np.pi**2/3./N  # k = {0, N/2}
        theory_var2 = np.pi**2/6./N  # otherwise
        for K in range(yk.size-1):
    #        aic[K] = np.sum(yk[K+1:]**2)/theory_var + 2.*K
            aic[K] = ( (1./theory_var2)*np.sum(yk[K+1:-1]**2) + (1./theory_var1)*yk[-1]**2 ) + 2.*(K+1)
        aic[-1] = 2.*yk.size
    else:
        aic[-1] = 0. # + (2*(yk.size+1))
        for K in range(yk.size-2, -1, -1):   # N-2, N-3, ..., 0
            aic[K] = aic[K+1] + yk[K+1]**2/theory_var[K+1]
        aic = aic + 2.*(np.arange(yk.size) + 1)
    return aic

def dct_AICc(yk, theory_var=None):
    """AICc[K] = AIC[K] + 2*(K+1)*(K+2)/(NF-K-2)
    Assumiamo di tenere tutti i k <= K."""
    aic = dct_AIC(yk, theory_var)
    KK = np.arange(yk.size-2)
    aic[:-2] = aic[:-2] + 2.*(KK+1)*(KK+2)/(yk.size-KK-2.)
    aic[-2] = aic[-3]  # not defined
    aic[-1] = aic[-3]  # not defined
    return aic

def dct_aic_ab(yk, theory_var, A=1.0, B=2.0):
    """AIC[K] = sum_{k>K} c_k^2/theory_var + 2*K
    Assumiamo di tenere tutti i k <= K."""
    aic = np.zeros(yk.size)
    aic[-1] = 0.
    for K in range(yk.size-2, -1, -1):   # N-2, N-3, ..., 0
        aic[K] = aic[K+1] + yk[K+1]**2/theory_var[K+1]
    aic = A*aic + B*(np.arange(yk.size) + 1)
    return aic

def dct_filter_psd(y, K=None):
    # K is the maximum coefficient summed (c_k = 0 per k > K)
    if (K >= y.size):
        print "! Warning:  dct_filter_psd K value ({:}) out of range.".format(K)
        return np.full(y.size, np.NaN)
    yk = dct(y, type=1)
    if K is not None:
        yk[K+1:] = 0.
    ynew = dct(yk, type=1)/(y.size-1)*0.5
    return ynew

def dct_filter_tau(y):
    # K is the maximum coefficient summed (c_k = 0 per k > K)
    yk = dct(y, type=1)/(y.size-1)
    ftau = np.zeros(y.size)
    ftau[0] = yk[0]
    ftau[1] = 0.5*yk[0] + yk[1]
    for i in range(2, yk.size-1):
        ftau[i] = ftau[i-1] + yk[i]
    ftau[-1] = ftau[-2] + 0.5*yk[-1]
    return ftau

#def cos_filter(y, K=None):
#    yk = dct(y, type=1) #*np.exp(-(xxx/50.)**8))#[:50]
#    if K is not None:
#        yk[K:] = 0.
#    ynew = dct(yk, type=1)/(y.size-1)*0.5
#    return ynew, yk
#
#def Cos_filter(y, K=None):      ## kept for compatibility purposes
#    yk = dct(y, type=1) #*np.exp(-(xxx/50.)**8))#[:50]
#    if K is not None:
#        yk[K:] = 0.
#    ynew = dct(yk, type=1)/(y.size-1)*0.5
#    return ynew, yk


#################################################################
#
# Produce_density takes as inputs the following numpy arrays:
#
#                INPUT :
#
# - aic: The vector with the estimate of the Akaike information aic(k), with k ranging from k_beg:k_max. 
#         I think it is important to have values till k_max but the first k_beg can be different from one ; 
# - sigma, mean : for the k_beg:k_max, in the same order, the mean and sigmas of the estimate
#                 of the transport coefficient ;        
# - method : a string equal to 'one' or 'two' with two different methods of inizializing p[k], the estimate probability
#            of a given k. I hope in real cases the two methods should be equivalent. 
#
#                RETURNS:
#
# The probability p[ik], a grid and a density for the probability of the transport coefficient on that grid, obtained as:
# density[igrid] \sim Sum_ik p[ik] * N[ mean[ik],sigma[ik] ].
# All pi factors are never inserted and taken care I hope by final normalization.
#
# The functions produce_p and produce_p_density distinguish the two conceptual steps of producing p and using this quantity 
# to provide the final estimate on the transport coefficient
#

def produce_p(aic, method='ba', force_normalize=False):
    k0 = np.argmin(aic)
    kM = len(aic)
    p  = np.zeros(kM)
    
    if (method == 'min'):
        # just the min(aic)
        p[k0] = 1.0
    
    elif (method == 'baroni'):
        for ik in range(kM):
            delta_aic = aic[ik] - aic[k0]  
            GAMMA = 3.9215536345675050924567623117545  # (epsilon/sigma)^2
            p[ik] = np.exp(-0.25 * np.log(1. + GAMMA) * delta_aic)

    elif ((method == 'burnham-anderson') or (method == 'ba')):
        delta_aic = aic - aic[k0]
        p = np.exp(-0.5 * delta_aic)
        p = p / np.sum(p)
    
    elif ((method == 'burnham-anderson2') or (method == 'ba2')):
        delta_aic = aic - aic[k0]
        p = np.exp(-delta_aic)
        p = p / np.sum(p)

    elif (method == 'two'):
        for ik in range(kM):
            delta_aic = aic[ik] - aic[k0]
            p[ik] = np.exp(-delta_aic ** 2 / ( 2.0 * ( kM - k0) )) / np.sqrt(kM - k0)

    elif (method == 'four'):
        for ik in range(kM):
            delta_aic = aic[ik] - aic[k0]
            p[ik] = np.exp(-delta_aic ** 2 / ( 2.0 * np.abs(ik - k0) ))
        p[k0] = 1.0

    else:
        raise KeyError('P_AIC METHOD not valid.')
    #p[ik] = np.exp(- delta_aic ** 2 / ( 2.0 * ( kM - ik) ) ) / np.sqrt(kM - ik)
    #p[ik] = np.exp(-delta_aic ** 2 / ( 2.0 * np.abs(ik - k0) )) / np.sqrt(np.abs(ik - k0))
    
    # normalize p
    if force_normalize:
        integral = np.trapz(p)
        for ik in range(kM):
            p[ik] = p[ik] / integral
    
    # checks
    if any(np.isnan(p)):
        raise Warning('WARNING: P contains NaN.')
    if (np.min(p < 0)):
        raise Warning('WARNING: P < 0.')
    return p


def produce_p_density(p, sigma, mean, grid=None, grid_size=1000):
    kM = len(mean)
    dm = np.min(mean)
    dM = np.max(mean)
    argdm = np.argmin(mean)
    argdM = np.argmax(mean)
    sigmam = sigma[argdm]
    sigmaM = sigma[argdM]
    if grid is None:
        return_grid = True
        # generate grid with grid_size points
        delta = ( (dM+5.*sigmaM) - (dm-5.*sigmam) ) / float(grid_size-1)
        grid = np.linspace( dm-5.*sigmam, dM+5.*sigmaM, grid_size )
    else:
        return_grid = False
        delta = grid[1] - grid[0]
    density = np.zeros(len(grid))
    for ik in range(kM):
        density = density + p[ik] * np.exp( -(grid-mean[ik])**2 / (2.0*(sigma[ik]**2)) ) / sigma[ik]
    somma = np.trapz(density) * delta
    density = density / somma
    if return_grid:
        return density, grid
    else:
        return density


def grid_statistics(grid, density, grid2=None):
    """Compute distribution mean and std.
      media   = \sum_i (density[i] * grid[i])
      std     = sqrt( \sum_i (density[i] * grid[i]^2) - media^2 )
       oppure = sqrt( \sum_i (density[i] * grid2[i])  - media^2 )"""
    somma = np.sum(density)
    media = np.dot(density, grid) / somma
    if grid2 is None:
        var = np.dot(density, grid**2) / somma - media**2
    else:
        var = np.dot(density, grid2) / somma - media**2
    std = np.sqrt(var)
    return media, std



################################################################################
################################################################################

class LowPassFilter(object):
    
    def __init__(self, *args, **kwargs):
        if (len(args) < 1):
            self.filtertype = kwargs.get('filtertype', None)
            if self.filtertype is None:
                raise ValueError('Not enough input arguments.')
        else:
            self.filtertype = args[0]
        if (self.filtertype == 'exp'):
            if (len(args) == 5):
                self.freqs  = args[1]
                self.f0     = args[2]
                self.alpha  = args[3]
                self.minatt = args[4]
            else:
                self.freqs  = kwargs.get('freqs', None)
                self.f0     = kwargs.get('f0', None)
                self.alpha  = kwargs.get('alpha', None)
                self.minatt = kwargs.get('minatt', None)
        return

        
    def __repr__(self):
        msg = '{} low-pass filter:\n'.format(self.filtertype) + \
              '  f0 = {}\n'.format(self.f0) + \
              '  alpha = {}\n'.format(self.alpha)
        return msg

    
    def compute_response(self, freqs=None):
        if freqs is not None:
            self.freqs = freqs
        if self.freqs is None:
            raise ValueError('freqs not set.')
        if (self.filtertype == 'exp'):
            self.response = self.exp_filter()
        self.logresponse = np.log(self.response)
        return
    
    
    def exp_filter(self):
        if self.f0 is None:
            raise ValueError('f0 not set.')
        if self.alpha is None:
            raise ValueError('alpha not set.')
        if self.minatt is None:
            self.minatt = 1.0e-3
        lf = np.zeros(self.freqs.size)
        for i in range(self.freqs.size):
            if (0. <= self.freqs[i] < 0.5):
                lf[i] = np.exp( -(self.freqs[i]/self.f0)**self.alpha )*(1.0 - self.minatt) + self.minatt
            elif (0.5 <= self.freqs[i] < 1.):
                lf[i] = np.exp( -(np.abs(self.freqs[i]-1.)/self.f0)**self.alpha )*(1.0 - self.minatt) + self.minatt
            else:
                print "ERROR: frequency out of range!"
        return lf

    
def resample_psd(freqs, psd, cutfrequency):
    if (cutfrequency >= freqs[-1]):
        return freqs, psd
    Nfreqs = freqs.size - 1
    cutidx = (np.abs(freqs - cutfrequency)).argmin()
    if (Nfreqs%cutidx == 0 ): # cut frequency is sub-multiple of max freq
        DT = Nfreqs/cutidx
        if (DT > 2):
            raise Warning('DT Not implemented.')
        newpsd = psd.copy()[:cutidx+1]
        newpsd = newpsd + psd[:-cutidx-2:-1]
        newfreqs = freqs[:cutidx+1]
        #print cutidx, DT, freqs[cutidx], newpsd.size
    else:
        raise Warning('Not implemented.')
    return newfreqs, newpsd
