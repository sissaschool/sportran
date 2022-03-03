# -*- coding: utf-8 -*-

import numpy as np


def freq_THz_to_red(f_THz, DT_FS):
    """
    Converts THz to reduced frequency units.
       f[red] = f[THz] * dt[fs] / 1000
    """
    return f_THz / 1000. * DT_FS


def freq_red_to_THz(f_red, DT_FS):
    """
    Converts reduced frequency units to THz.
       f[THz] = f[red] * 1000 / dt[fs]
    """
    return f_red * 1000. / DT_FS


def logtau_to_tau(logtau, logtau_mean, logtau_var, correct_mean=True):
    if correct_mean:
        logtau = logtau - logtau_mean
        tau = np.exp(logtau)
        tau_std = np.sqrt(tau * logtau_var)
    else:
        raise NotImplementedError()
    return tau, tau_std


def generate_empirical_spectrum(psd):
    """Add noise to a periodogram and generate a complex spectrum."""
    ### Probabilmente ci vuole un fattore 1/N, in maniera da tirar via il fattore N
    ### in compute_trajectory, cosi' diventa consistente con compute_spectrum
    spectr = np.random.normal(0., np.sqrt(0.5*psd*(psd.size-1)*2)) + \
               1.0J*np.random.normal(0., np.sqrt(0.5*psd*(psd.size-1)*2))
    spectr[0] = spectr[0].real
    spectr[-1] = spectr[-1].real
    return spectr
