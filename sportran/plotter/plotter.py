# -*- coding: utf-8 -*-
"""
Defines an (abstract) Plotter class and all the plot functions that its subclasses can import.
"""

import os
import math
import numpy as np
from sportran.utils import log
from sportran.md.tools.spectrum import freq_THz_to_red

# import the matplotlib pyplot module loaded by the __init__
from . import plt
from matplotlib.ticker import MultipleLocator
from matplotlib.backends.backend_pdf import PdfPages

# list of colors
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

################################################################################


class Plotter():
    """
    Plotter abstract class. Essentially empty.
    """
    _plot_style = None
    pass


################################################################################


def _n_tick_in_range(beg, end, n):
    size = end - beg
    n_cifre = math.floor(math.log(size / n, 10.0))
    delta = math.ceil((size / n) / 10**n_cifre) * 10**n_cifre
    return delta, delta / 2


def _index_cumsum(arr, p):
    if p > 1 or p < 0:
        raise ValueError('p must be between 0 and 1')
    arr_int = np.cumsum(arr)
    arr_int = arr_int / arr_int[-1]
    idx = 0
    while arr_int[idx] < p:
        idx = idx + 1
    return idx


def addPlotToPdf(func, pdf, *args, **kwargs):
    result = func(*args, **kwargs)
    pdf.savefig()
    plt.close()
    return result


################################################################################
## Plot functions
## The first argument should be a Current or MDSample object, if they are supposed to be transformed into a method by the add_method decorator.


def plot_trajectory(x, *, axes=None, FIGSIZE=None, **plot_kwargs):
    """
    Plot the time series.
    """
    if x.traj is None:
        raise ValueError('Trajectory not defined.')
    if axes is None:
        figure, axes = plt.subplots(1, figsize=FIGSIZE)
    axes.plot(x.traj, **plot_kwargs)
    axes.set_xlabel(r'$t$ [ps]')
    axes.grid()
    return axes

def plot_periodogram(current, PSD_FILTER_W=None, *, freq_units='THz', freq_scale=1.0, axes=None, kappa_units=False,
                     FIGSIZE=None, mode='log', **plot_kwargs):   # yapf: disable
    """
    Plots the current's periodogram (psd)
    :param current:         current object to plot periodogram
    :param PSD_FILTER_W:    width of the filtering window
    :param freq_units:      'thz'  [THz]
                            'red'  [omega*DT/(2*pi)]
    :param freq_scale:      rescale red frequencies by this factor (e.g. 2 --> freq = [0, 0.25])
    :param axes:            plot periodograms in units of kappa (default: False) - NB: log-psd not converted
    :param kappa_units:     matplotlib.axes.Axes object (if None, create one)
    :param FIGSIZE:         size of the plot

    :return: a matplotlib.axes.Axes object
    """

    # recompute PSD if needed
    if current.psd is None:
        current.compute_psd()
    # (re)compute filtered psd, if a window has been defined
    if (PSD_FILTER_W is not None) or (current.PSD_FILTER_W is not None):
        current.filter_psd(PSD_FILTER_W, freq_units)
    else:   # use a zero-width (non-filtering) window
        current.filter_psd(0.)
    if kappa_units:   # plot psd in units of kappa - the log(psd) is not converted
        psd_scale = 0.5 * current.KAPPA_SCALE
    else:
        psd_scale = 1.0

    if axes is None:
        figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
        plt.subplots_adjust(hspace=0.1)
    if freq_units in ('THz', 'thz'):
        axes[0].plot(current.freqs_THz, psd_scale * current.fpsd, **plot_kwargs)
        axes[0].set_xlim([0., current.Nyquist_f_THz])
        if mode == 'log':
            axes[1].plot(current.freqs_THz, current.flogpsd, **plot_kwargs)
            axes[1].set_xlim([0., current.Nyquist_f_THz])
            axes[1].set_xlabel(r'$f$ [THz]')
    elif freq_units == 'red':
        axes[0].plot(current.freqs / freq_scale, psd_scale * current.fpsd, **plot_kwargs)
        axes[0].set_xlim([0., 0.5 / freq_scale])
        if mode == 'log':
            axes[1].plot(current.freqs / freq_scale, current.flogpsd, **plot_kwargs)
            axes[1].set_xlim([0., 0.5 / freq_scale])
            axes[1].set_xlabel(r'$f$ [$\omega$*DT/2$\pi$]')
    else:
        raise ValueError('Frequency units not valid.')
    axes[0].xaxis.set_ticks_position('top')
    if kappa_units:
        axes[0].set_ylabel(r'PSD [{}]'.format(current._KAPPA_SI_UNITS))
    else:
        axes[0].set_ylabel(r'PSD')
    axes[0].grid()
    if mode == 'log':
        axes[1].xaxis.set_ticks_position('bottom')
        axes[1].set_ylabel(r'log(PSD)')
        axes[1].grid()
    #print('axes2=', axes)
    return axes


def plot_ck(current, *, axes=None, label=None, FIGSIZE=None):
    """
    Plots the cepstral coefficients c_K.
    :param current: current object to plot
    :param axes: matplotlib.axes.Axes object (if None, create one)
    :param label:
    :param FIGSIZE: size of the plot

    :return: a matplotlib.axes.Axes object
    """

    if axes is None:
        figure, axes = plt.subplots(1, figsize=FIGSIZE)
    color = next(axes._get_lines.prop_cycler)['color']
    axes.plot(current.dct.logpsdK, 'o-', c=color, label=label)
    axes.plot(current.dct.logpsdK + current.dct.logpsdK_THEORY_std, '--', c=color)
    axes.plot(current.dct.logpsdK - current.dct.logpsdK_THEORY_std, '--', c=color)
    axes.axvline(x=current.dct.aic_Kmin, ls='--', c=color)
    axes.set_xlabel(r'$k$')
    axes.set_ylabel(r'$c_k$')
    return axes


def plot_L0_Pstar(current, *, axes=None, label=None, FIGSIZE=None):
    """
    Plots L0 as a function of P*.
    :param current:         current object to plot
    :param axes:            matplotlib.axes.Axes object (if None, create one)
    :param label:
    :param FIGSIZE:         size of the plot

    :return: a matplotlib.axes.Axes object
    """
    if axes is None:
        figure, axes = plt.subplots(1, figsize=FIGSIZE)
    color = next(axes._get_lines.prop_cycler)['color']
    axes.plot(np.arange(current.NFREQS) + 1, current.dct.logtau, '.-', c=color, label=label)
    axes.plot(np.arange(current.NFREQS) + 1, current.dct.logtau + current.dct.logtau_THEORY_std, '--', c=color)
    axes.plot(np.arange(current.NFREQS) + 1, current.dct.logtau - current.dct.logtau_THEORY_std, '--', c=color)
    axes.axvline(x=current.dct.aic_Kmin + 1, ls='--', c=color)
    axes.set_xlim([0, 3 * current.dct.aic_Kmin])
    max_y = np.amax((current.dct.logtau + current.dct.logtau_THEORY_std)[current.dct.aic_Kmin:3 * current.dct.aic_Kmin])
    min_y = np.amin((current.dct.logtau - current.dct.logtau_THEORY_std)[current.dct.aic_Kmin:3 * current.dct.aic_Kmin])
    axes.set_ylim([min_y * 0.8, max_y * 1.2])
    axes.set_xlabel(r'$P^*$')
    axes.set_ylabel(r'$L_0(P*)$')
    return axes


def plot_kappa_Pstar(current, *, axes=None, label=None, FIGSIZE=None):
    """
    Plots the value of kappa as a function of P*.
    :param current: current object to plot
    :param axes: matplotlib.axes.Axes object (if None, create one)
    :param label:
    :param FIGSIZE: size of the plot

    :return: a matplotlib.axes.Axes object
    """
    if axes is None:
        figure, axes = plt.subplots(1, figsize=FIGSIZE)
    color = next(axes._get_lines.prop_cycler)['color']
    axes.plot(np.arange(current.NFREQS) + 1, current.dct.tau * current.KAPPA_SCALE * 0.5, '.-', c=color, label=label)
    axes.plot(np.arange(current.NFREQS) + 1, (current.dct.tau + current.dct.tau_THEORY_std) * current.KAPPA_SCALE * 0.5,
              '--', c=color)   # yapf: disable
    axes.plot(np.arange(current.NFREQS) + 1, (current.dct.tau - current.dct.tau_THEORY_std) * current.KAPPA_SCALE * 0.5,
              '--', c=color)   # yapf: disable
    axes.axvline(x=current.dct.aic_Kmin + 1, ls='--', c=color)
    axes.axhline(y=current.kappa_Kmin, ls='--', c=color)
    axes.set_xlim([0, 3 * current.dct.aic_Kmin])
    max_y = np.amax(current.KAPPA_SCALE * 0.5 *
                    (current.dct.tau + current.dct.tau_THEORY_std)[current.dct.aic_Kmin:3 * current.dct.aic_Kmin])
    min_y = np.amin(current.KAPPA_SCALE * 0.5 *
                    (current.dct.tau - current.dct.tau_THEORY_std)[current.dct.aic_Kmin:3 * current.dct.aic_Kmin])
    axes.set_ylim([min_y * 0.8, max_y * 1.2])
    axes.set_xlabel(r'$P^*$')
    axes.set_ylabel(r'$\kappa(P^*)$ [{}]'.format(current._KAPPA_SI_UNITS))
    return axes

def plot_cepstral_spectrum(current, *, freq_units='THz', freq_scale=1.0, axes=None, kappa_units=True, FIGSIZE=None, mode='log',
                           **plot_kwargs):   # yapf: disable
    """
    Plots the cepstral spectrum.

    :param current:         current object to plot
    :param freq_units:      'thz'  [THz]
                            'red'  [omega*DT/(2*pi)]
    :param freq_scale:      rescale red frequencies by this factor (e.g. 2 --> freq = [0, 0.25])
    :param axes:            plot periodograms in units of kappa (default: False) - NB: log-psd not converted
    :param kappa_units:     matplotlib.axes.Axes object (if None, create one)
    :param FIGSIZE:         size of the plot

    :return: a matplotlib.axes.Axes object
    """
    if axes is None:
        figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
    plt.subplots_adjust(hspace=0.1)
    if kappa_units:
        psd_scale = 0.5 * current.KAPPA_SCALE
    else:
        psd_scale = 1.0
    if freq_units in ('THz', 'thz'):
        axes[0].plot(current.freqs_THz, current.dct.psd * psd_scale, **plot_kwargs)
        axes[0].set_xlim([0., current.Nyquist_f_THz])
        if mode == 'log':
            axes[1].plot(current.freqs_THz, current.dct.logpsd, **plot_kwargs)
            axes[1].set_xlim([0., current.Nyquist_f_THz])
            axes[1].set_xlabel(r'$f$ [THz]')
    elif freq_units == 'red':
        axes[0].plot(current.freqs / freq_scale, current.dct.psd * psd_scale, **plot_kwargs)
        axes[0].set_xlim([0., 0.5 / freq_scale])
        if mode == 'log':
            axes[1].plot(current.freqs / freq_scale, current.dct.logpsd, **plot_kwargs)
            axes[1].set_xlim([0., 0.5 / freq_scale])
            axes[1].set_xlabel(r'$f$ [$\omega$*DT/2$\pi$]')
    else:
        raise ValueError('Units not valid.')
    axes[0].xaxis.set_ticks_position('top')
    axes[0].set_ylabel(r'PSD')
    if kappa_units:
        axes[0].set_ylabel(r'PSD [{}]'.format(current._KAPPA_SI_UNITS))
    else:
        axes[0].set_ylabel(r'PSD')
    axes[0].grid()
    if mode == 'log':
        axes[1].xaxis.set_ticks_position('bottom')
        axes[1].set_ylabel(r'log(PSD)')
        axes[1].grid()
    return axes


def plot_fstar_analysis(current, xf, FSTAR_THZ_LIST, *, axes=None, FIGSIZE=None, **plot_kwargs):
    """
    Plots kappa(P*) as a function of the f*.
    """
    if axes is None:
        figure, ax = plt.subplots(2, sharex=True, figsize=FIGSIZE)
    else:
        ax = axes
    ax[0].errorbar(FSTAR_THZ_LIST, [xff.kappa_Kmin for xff in xf], yerr=[xff.kappa_Kmin_std for xff in xf],
                   **plot_kwargs)
    ax[1].errorbar(FSTAR_THZ_LIST, [xff.dct.logtau_Kmin for xff in xf], yerr=[xff.dct.logtau_std_Kmin for xff in xf],
                   **plot_kwargs)
    # ax[0].plot(current.freqs_THz, current.fpsd,    **plot_kwargs)
    # ax[1].plot(current.freqs_THz, current.flogpsd, **plot_kwargs)
    ax[0].xaxis.set_ticks_position('top')
    ax[0].set_ylabel(r'PSD')
    ax[0].grid()
    ax[1].xaxis.set_ticks_position('bottom')
    ax[1].set_xlabel(r'$f$ [THz]')
    ax[1].set_ylabel(r'log(PSD)')
    ax[1].grid()

    if axes is None:
        ax2 = [ax[0].twinx(), ax[1].twinx()]
        color = next(ax[0]._get_lines.prop_cycler)['color']
        color = next(ax[1]._get_lines.prop_cycler)['color']
        plot_periodogram(current, axes=ax2, c=color)
        ax[0].set_ylabel(r'$\kappa$ [W/(m*K)]')
        ax[1].set_ylabel(r'$\kappa$ [W/(m*K)]')
        return xf, ax, figure
    else:
        return xf, ax


def plot_resample(x, xf, PSD_FILTER_W=None, *, freq_units='THz', axes=None, FIGSIZE=None, mode='log'):
    """
    Plots the periodogram of a time series and of a filtered/resampled one for comparison.
    :param x:               a time series object to plot
    :param xf:              a filtered & resampled time series object
    :param freq_units:      'thz'  [THz]
                            'red'  [omega*DT/(2*pi)]
    :param PSD_FILTER_W:    PSD filtering window width [chosen frequency units]
    :param FIGSIZE:         plot figure size

    :return:                xf: a filtered & resampled time series object
                             axes: a matplotlib.axes.Axes object

    """
    fstar_THz = xf.Nyquist_f_THz
    TSKIP = int(x.Nyquist_f_THz / xf.Nyquist_f_THz)

    if not axes:
        figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
        axes = plot_periodogram(x, PSD_FILTER_W=PSD_FILTER_W, freq_units=freq_units, axes=axes, mode=mode,
                                kappa_units=True)   # this also updates x.PSD_FILTER_W
    xf.plot_periodogram(freq_units=freq_units, freq_scale=TSKIP, axes=axes, mode=mode, kappa_units=True)
    if freq_units in ('THz', 'thz'):
        axes[0].axvline(x=fstar_THz, ls='--', c='k')
        axes[0].set_xlim([0., x.Nyquist_f_THz])
        if mode == 'log':
            axes[1].axvline(x=fstar_THz, ls='--', c='k')
            axes[1].set_xlim([0., x.Nyquist_f_THz])
    elif freq_units == 'red':
        axes[0].axvline(x=0.5 / TSKIP, ls='--', c='k')
        axes[0].set_xlim([0., 0.5])
        if mode == 'log':
            axes[1].axvline(x=0.5 / TSKIP, ls='--', c='k')
            axes[1].set_xlim([0., 0.5])
    return axes


################################################################################
## DUPLICATE FUNCTIONS THAT NEED TO BE MERGED IF POSSIBLE


def plot_psd(jf, j2=None, j2pl=None, f_THz_max=None, k_SI_max=None, k_tick=None, f_tick=None):
    if f_THz_max is None:
        idx_max = _index_cumsum(jf.psd, 0.95)
        f_THz_max = jf.freqs_THz[idx_max]
    else:
        maxT = jf.freqs_THz[-1]
        if j2 is not None:
            if j2.freqs_THz[-1] > maxT:
                maxT = j2.freqs_THz[-1]
        if j2pl is not None:
            if j2pl.freqs_THz[-1] > maxT:
                maxT = j2pl.freqs_THz[-1]
        if maxT < f_THz_max:
            f_THz_max = maxT

    if k_SI_max is None:
        k_SI_max = np.max(
            jf.fpsd[:int(jf.freqs_THz.shape[0] * f_THz_max / jf.freqs_THz[-1])] * jf.KAPPA_SCALE * 0.5) * 1.3

    figure, ax = plt.subplots(1, 1)   # figsize=(3.8, 2.3)
    ax.plot(jf.freqs_THz, jf.psd * jf.KAPPA_SCALE * 0.5, lw=0.2, c='0.8', zorder=0)
    ax.plot(jf.freqs_THz, jf.fpsd * jf.KAPPA_SCALE * 0.5, c=colors[0], zorder=2)
    if j2 is not None:
        plt.axvline(x=j2.Nyquist_f_THz, ls='--', c='k', dashes=(1.4, 0.6), zorder=3)
    if j2pl is not None:
        plt.plot(j2pl.freqs_THz, j2pl.dct.psd * j2pl.KAPPA_SCALE * 0.5, c=colors[1], zorder=1)
    try:
        plt.plot(jf.freqs_THz, np.real(jf.fcospectrum[0][0]) * jf.KAPPA_SCALE * 0.5, c=colors[3], lw=1.0, zorder=1)
    except:
        pass

    ax.set_ylim([0, k_SI_max])
    ax.set_xlim([0, f_THz_max])
    ax.set_xlabel(r'$\omega/2\pi$ (THz)')
    ax.set_ylabel(r'${{}}^{{\ell}}\hat{{S}}_{{\,k}}$ [{}]'.format(jf._KAPPA_SI_UNITS))

    if f_tick is None:
        dx1, dx2 = _n_tick_in_range(0, f_THz_max, 5)
    else:
        dx1 = f_tick
        dx2 = dx1 / 2
    if k_tick is None:
        dy1, dy2 = _n_tick_in_range(0, k_SI_max, 5)
    else:
        dy1 = k_tick
        dy2 = dy1 / 2

    ax.xaxis.set_major_locator(MultipleLocator(dx1))
    ax.xaxis.set_minor_locator(MultipleLocator(dx2))
    ax.yaxis.set_major_locator(MultipleLocator(dy1))
    ax.yaxis.set_minor_locator(MultipleLocator(dy2))


################################################################################
## MAYBE OBSOLETE STUFF


def plot_cepstral_conv(jf, pstar_max=None, k_SI_max=None, pstar_tick=None, kappa_tick=None):
    """
    MAYBE OBSOLETE
    """
    if pstar_max is None:
        pstar_max = (jf.dct.aic_Kmin + 1) * 2.5
    if k_SI_max is None:
        k_SI_max = jf.dct.tau[jf.dct.aic_Kmin] * jf.KAPPA_SCALE

    f, ax2 = plt.subplots(1, 1)   # figsize=(3.8, 2.3)
    ax2.axvline(x=jf.dct.aic_Kmin + 1, ls='--', c='k', dashes=(1.4, 0.6), zorder=-3)
    ax2.fill_between(
        np.arange(jf.dct.logtau.shape[0]) + 1, (jf.dct.tau - jf.dct.tau_THEORY_std) * jf.KAPPA_SCALE * 0.5,
        (jf.dct.tau + jf.dct.tau_THEORY_std) * jf.KAPPA_SCALE * 0.5, alpha=0.3, color=colors[4], zorder=-3)
    ax2.plot(
        np.arange(jf.dct.logtau.shape[0]) + 1, jf.dct.tau * jf.KAPPA_SCALE * 0.5, label=r'Cepstral method', marker='o',
        c=colors[4], zorder=-3)
    ax2.set_xlabel(r'$P^*$')
    ax2.set_ylabel(r'$\kappa$ [{}]'.format(jf._KAPPA_SI_UNITS))
    ax2.set_xlim([0, pstar_max])
    ax2.set_ylim([0, k_SI_max])
    # ax2.grid()
    ax2.legend()

    if pstar_tick is None:
        dx1, dx2 = _n_tick_in_range(0, pstar_max, 5)
    else:
        dx1 = pstar_tick
        dx2 = dx1 / 2
    if kappa_tick is None:
        dy1, dy2 = _n_tick_in_range(0, k_SI_max, 5)
    else:
        dy1 = kappa_tick
        dy2 = dy1 / 2
    ax2.xaxis.set_major_locator(MultipleLocator(dx1))
    ax2.xaxis.set_minor_locator(MultipleLocator(dx2))
    ax2.yaxis.set_major_locator(MultipleLocator(dy1))
    ax2.yaxis.set_minor_locator(MultipleLocator(dy2))


def plot_other(jf, idx1, idx2, f_THz_max=None, k_SI_max=None, k_SI_min=None, k_tick=None, f_tick=None):
    """
    MAYBE OBSOLETE
    """
    if f_THz_max is None:
        idx_max = _index_cumsum(np.abs(jf.fcospectrum[idx1][idx2]), 0.95)
        f_THz_max = jf.freqs_THz[idx_max]
    else:
        maxT = jf.freqs_THz[-1]
        if maxT < f_THz_max:
            f_THz_max = maxT

    if k_SI_max is None:
        k_SI_max = np.max(
            np.abs(jf.fcospectrum[idx1][idx2])[:int(jf.freqs_THz.shape[0] * f_THz_max / jf.freqs_THz[-1])] *
            jf.KAPPA_SCALE * 0.5) * 1.3
    if k_SI_min is None:
        k_SI_min = -k_SI_max

    figure, ax = plt.subplots(1, 1)   # figsize=(3.8, 2.3)
    ax.plot(jf.freqs_THz, np.real(jf.fcospectrum[idx1][idx2]) * jf.KAPPA_SCALE * 0.5, c=colors[3], lw=1.0, zorder=1)
    ax.plot(jf.freqs_THz, np.imag(jf.fcospectrum[idx1][idx2]) * jf.KAPPA_SCALE * 0.5, c=colors[2], lw=1.0, zorder=1)

    ax.set_ylim([k_SI_min, k_SI_max])
    ax.set_xlim([0, f_THz_max])
    ax.set_xlabel(r'$\omega/2\pi$ (THz)')
    ax.set_ylabel(r'$S^{{{}{}}}$'.format(idx1, idx2))

    if f_tick is None:
        dx1, dx2 = _n_tick_in_range(0, f_THz_max, 5)
    else:
        dx1 = f_tick
        dx2 = dx1 / 2
    if k_tick is None:
        dy1, dy2 = _n_tick_in_range(0, k_SI_max, 5)
    else:
        dy1 = k_tick
        dy2 = dy1 / 2

    ax.xaxis.set_major_locator(MultipleLocator(dx1))
    ax.xaxis.set_minor_locator(MultipleLocator(dx2))
    ax.yaxis.set_major_locator(MultipleLocator(dy1))
    ax.yaxis.set_minor_locator(MultipleLocator(dy2))
