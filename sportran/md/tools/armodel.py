# -*- coding: utf-8 -*-

import numpy as np


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
            self.P = kwargs.get('order', None)
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
        self.phi_std = np.sqrt(np.diag(self.cov)[:-1])
        self.sigma2_std = np.sqrt(self.cov[-1, -1])
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
                x_t = x_t + traj[t - i - 1] * self.phi[i]
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
    a = np.zeros((P, P))
    for i in range(P):
        for j in range(i, P):
            c = 0.
            for t in range(P, RUN_TIME):
                c += y[t - (i + 1)] * y[t - (j + 1)] / (RUN_TIME - P)
            a[i, j] = a[j, i] = c

    b = np.zeros(P)
    for i in range(P):
        c = 0.
        for t in range(P, RUN_TIME):
            c += y[t - (i + 1)] * y[t] / (RUN_TIME - P)
        b[i] = c

    phi = solve(a, b, sym_pos=True)

    # compute sigma2
    sigma2 = 0.
    for t in range(P, RUN_TIME):
        sigma2 += (y[t] - np.dot(y[t - P:t], phi[::-1]))**2
        #res = y[t]
        #for j in range(P):
        #    res -= phi[j]*y[t-(j+1)]
        #sigma2 += res**2
    sigma2 *= 1. / (RUN_TIME - P)

    # compute asymptotic covariance matrix of parameters
    yy = y[P - 1:-1]
    for i in range(P - 1):
        yy = np.row_stack((yy, y[P - 2 - i:-2 - i]))

    V = np.zeros((P + 1, P + 1))
    V[:P, :P] = inv(np.cov(yy)) * sigma2 / (RUN_TIME - P - 1)
    V[P, P] = 2. * sigma2 * sigma2 / (RUN_TIME - P - 1)
    return phi, sigma2, V


def ar_psd(AR_phi, AR_sigma2, NFREQS, DT=1):
    """Compute psd of an AR(P) process. BE CAREFUL WITH NORMALIZATION IF DT!=1"""
    P = len(AR_phi)
    freqs = np.linspace(0., 0.5 / DT, NFREQS + 1)
    AR_psd = np.zeros(freqs.size)
    for i in range(freqs.size):
        phiz = np.sum(AR_phi * np.exp(-2.0j * np.pi * freqs[i] * np.arange(1, P + 1)))
        #AR_psd[i] = 2. * DT * AR_sigma2 / np.abs( 1. - phiz )**2  # factor 2 comes from one-sided
        AR_psd[i] = DT * AR_sigma2 / np.abs(1. - phiz)**2
    return freqs / DT, AR_psd


def ar_tau(AR_phi, AR_sigma2, AR_phi_cov=None, DT=1, RUN_TIME=0):
    """Compute tau of an AR(P) process."""

    P = AR_phi.size

    if AR_phi_cov is not None:
        # check if covariance matrix contains cov(sigma2,_)
        if (AR_phi_cov.shape == (P + 1, P + 1)):
            AR_cov = AR_phi_cov
        elif AR_phi_cov.shape == (P, P):
            if (RUN_TIME == 0):
                raise ValueError('You should pass RUN_TIME to AR_tau function.')
            # build total covariance matrix (phi+sigma2)
            AR_cov = np.zeros((P + 1, P + 1))
            AR_cov[:-1, :-1] = AR_phi_cov
            AR_cov[-1, -1] = 2. * AR_sigma2**2 / (RUN_TIME - P - 1)
        else:
            raise ValueError('AR_phi_cov matrix has wrong dimensions.')

    # compute model tau from S(f=0)
    phiz = np.sum(AR_phi)
    AR_tau = 0.5 * DT * AR_sigma2 / np.abs(1. - phiz)**2

    if AR_phi_cov is not None:
        grad_tau = DT * np.append(AR_sigma2 / phiz**3 * np.ones(P), 0.5 / phiz**2)
        AR_tau_stderr = np.sqrt(grad_tau.dot(AR_cov).dot(grad_tau))
        return AR_tau, AR_tau_stderr
    else:
        return AR_tau
