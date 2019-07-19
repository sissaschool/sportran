import numpy as np
from scipy.signal import lfilter
from thermocepstrum.utils.utils import PrintMethod
log = PrintMethod()

################################################################################


def integrate_acf(acf):
    """Returns the integral function of acf, i.e. its integral as a function of the upper integration limit.
    Supports multi-component (N, N_COMP) arrays.
    Trapezial integration is used.
          tau[i] = trapz_{0}^{i} acf
    """

    N = acf.shape[0]
    tau = np.zeros(acf.shape)
    for i in range(1, N):
        tau[i] = tau[i - 1] + 0.5 * acf[i - 1] + 0.5 * acf[i]
    return tau


def runavefilter(X, WF):
    """Computes the running average of a numpy array over WF consecutive elements (or WF+1 if WF is even):
            (X[i-WF/2]+...+X[i]+...+X[i+WF/2]) / WF
    assumes that the array is "even" ( X[-i] = X[i] ) and anti-periodic (X[(N-1)+i]=X[(N-1)-i]), like a ONE-SIDED PSD.
    """

    if (WF % 2 == 0):
        WF = WF + 1
    W = int(WF / 2)

    Y = np.concatenate((X[W:0:-1], X, X[-2:-W - 2:-1]))
    return np.convolve(Y, np.array([1.0 / WF] * WF), 'valid')


################################################################################


def filter_and_sample(y_big, W, DT, window='rectangular', even_NSTEPS=True, detrend=False, drop_first=True):
    """Filter signal with moving average window of width W and then sample it
    with time step DT."""

    if (W > 1):
        if (window == 'rectangular'):
            y_f = lfilter((1. / W) * np.ones(W), 1., y_big, axis=0)
        else:
            raise NotImplementedError('Not implemented window type.')

        # drop first W steps (initial conditions)
        if drop_first:
            y = y_f[(W - 1)::DT]
        else:
            y = y_f[::DT]
    else:
        y = y_big[::DT]

    # remove the mean
    if detrend:
        y = y - np.mean(y, axis=0)

    # keeps an even number of points
    if even_NSTEPS:
        if (y.shape[0] % 2 == 1):
            return y[:-1]
    return y


def generate_empirical_spectrum(psd):
    """Add noise to a periodogram and generate a complex spectrum."""
    ### Probabilmente ci vuole un fattore 1/N, in maniera da tirar via il fattore N
    ### in compute_trajectory, cosi' diventa consistente con compute_spectrum
    spectr = np.random.normal(0., np.sqrt(0.5*psd*(psd.size-1)*2)) + \
               1.0J*np.random.normal(0., np.sqrt(0.5*psd*(psd.size-1)*2))
    spectr[0] = spectr[0].real
    spectr[-1] = spectr[-1].real
    return spectr


################################################################################


def logtau_to_tau(logtau, logtau_mean, logtau_var, correct_mean=True):
    if correct_mean:
        logtau = logtau - logtau_mean
        tau = np.exp(logtau)
        tau_std = np.sqrt(tau * logtau_var)
    else:
        raise NotImplementedError()
    return tau, tau_std


def resample_psd(freqs, psd, cutfrequency):
    if (cutfrequency >= freqs[-1]):
        return freqs, psd
    Nfreqs = freqs.size - 1
    cutidx = (np.abs(freqs - cutfrequency)).argmin()
    if (Nfreqs % cutidx == 0):   # cut frequency is sub-multiple of max freq
        DT = Nfreqs / cutidx
        if (DT > 2):
            raise Warning('DT Not implemented.')
        newpsd = psd.copy()[:cutidx + 1]
        newpsd = newpsd + psd[:-cutidx - 2:-1]
        newfreqs = freqs[:cutidx + 1]
        #log.write_log(cutidx, DT, freqs[cutidx], newpsd.size)
    else:
        raise NotImplementedError('Not implemented.')
    return newfreqs, newpsd
