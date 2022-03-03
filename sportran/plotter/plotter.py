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


def plot_trajectory(x, *, axis=None, FIGSIZE=None, **plot_kwargs):
    """
    Plot the time series.
    """
    if x.traj is None:
        raise ValueError('Trajectory not defined.')
    if axis is None:
        figure, axis = plt.subplots(1, figsize=FIGSIZE)
    axis.plot(x.traj, **plot_kwargs)
    axis.set_xlabel(r'$t$ [ps]')
    axis.grid()
    return axis


def plot_periodogram(current, PSD_FILTER_W=None, *, freq_units='THz', freq_scale=1.0, axes=None, kappa_units=True,
                     FIGSIZE=None, mode='log', **plot_kwargs):
    """
    Plots the current's periodogram (psd)
    :param current:         current object to plot periodogram
    :param PSD_FILTER_W:    width of the filtering window
    :param freq_units:      'thz'  [THz]
                            'red'  [omega*DT/(2*pi)]
    :param freq_scale:      rescale red frequencies by this factor (e.g. 2 --> freq = [0, 0.25])
    :param axes:            matplotlib.axes.Axes object (if None, create one)
    :param kappa_units:     plot periodograms in units of kappa (default: True) - NB: log-psd not converted
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
    return axes


def plot_cospectrum_component(current, idx1, idx2, *, axis=None, FIGSIZE=None, f_THz_max=None, k_SI_max=None,
                              k_SI_min=None, k_tick=None, f_tick=None):
    """
    Plot the (idx1, idx2) component of the cospectrum.
    """
    if axis is None:
        figure, axis = plt.subplots(1, figsize=FIGSIZE)
    color1 = next(axis._get_lines.prop_cycler)['color']
    color2 = next(axis._get_lines.prop_cycler)['color']
    axis.plot(current.freqs_THz, np.real(current.fcospectrum[idx1][idx2]) * current.KAPPA_SCALE * 0.5, c=color1)
    axis.plot(current.freqs_THz, np.imag(current.fcospectrum[idx1][idx2]) * current.KAPPA_SCALE * 0.5, c=color2)

    if f_THz_max is None:
        f_THz_max = current.freqs_THz[_index_cumsum(np.abs(current.fcospectrum[idx1][idx2]), 0.95)]
    else:
        f_THz_max = min(f_THz_max, current.freqs_THz[-1])
    axis.set_xlim([0, f_THz_max])
    if k_SI_max is None:
        k_SI_max = np.max(
            np.abs(current.fcospectrum[idx1][idx2])[:int(current.NFREQS * f_THz_max / current.freqs_THz[-1])] *
            current.KAPPA_SCALE * 0.5) * 1.3
    if k_SI_min is None:
        k_SI_min = -k_SI_max
    axis.set_ylim([k_SI_min, k_SI_max])
    axis.set_xlabel(r'$\omega/2\pi$ (THz)')
    axis.set_ylabel(r'$S^{{{}{}}}$'.format(idx1, idx2))

    if f_tick is None:
        dx1, dx2 = _n_tick_in_range(0, f_THz_max, 5)
    else:
        dx1, dx2 = (f_tick, f_tick / 2)
    if k_tick is None:
        dy1, dy2 = _n_tick_in_range(0, k_SI_max, 5)
    else:
        dy1, dy2 = (k_tick, k_tick / 2)

    axis.xaxis.set_major_locator(MultipleLocator(dx1))
    axis.xaxis.set_minor_locator(MultipleLocator(dx2))
    axis.yaxis.set_major_locator(MultipleLocator(dy1))
    axis.yaxis.set_minor_locator(MultipleLocator(dy2))


def plot_ck(current, *, axis=None, label=None, FIGSIZE=None):
    """
    Plots the cepstral coefficients c_K.
    :param current: current object to plot
    :param axis: matplotlib.axes.Axes object (if None, create one)
    :param label:
    :param FIGSIZE: size of the plot

    :return: a matplotlib.axis.Axes object
    """

    if axis is None:
        figure, axis = plt.subplots(1, figsize=FIGSIZE)
    color = next(axis._get_lines.prop_cycler)['color']
    axis.plot(current.cepf.logpsdK, 'o-', c=color, label=label)

    axis.plot(current.cepf.logpsdK + current.cepf.logpsdK_THEORY_std, '--', c=color)
    axis.plot(current.cepf.logpsdK - current.cepf.logpsdK_THEORY_std, '--', c=color)
    axis.axvline(x=current.cepf.aic_Kmin, ls=':', c=color)
    axis.axvline(x=current.cepf.cutoffK, ls='--', c=color)
    axis.set_xlabel(r'$k$')
    axis.set_ylabel(r'$c_k$')
    return axis


def plot_L0_Pstar(current, *, axis=None, label=None, FIGSIZE=None):
    """
    Plots L0 as a function of P*.
    :param current:         current object to plot
    :param axis:            matplotlib.axes.Axes object (if None, create one)
    :param label:
    :param FIGSIZE:         size of the plot

    :return: a matplotlib.axis.Axes object
    """
    if axis is None:
        figure, axis = plt.subplots(1, figsize=FIGSIZE)
    color = next(axis._get_lines.prop_cycler)['color']
    axis.plot(np.arange(current.NFREQS) + 1, current.cepf.logtau, '.-', c=color, label=label)
    axis.plot(np.arange(current.NFREQS) + 1, current.cepf.logtau + current.cepf.logtau_THEORY_std, '--', c=color)
    axis.plot(np.arange(current.NFREQS) + 1, current.cepf.logtau - current.cepf.logtau_THEORY_std, '--', c=color)
    axis.axvline(x=current.cepf.aic_Kmin + 1, ls=':', c=color)
    axis.axvline(x=current.cepf.cutoffK + 1, ls='--', c=color)
    axis.set_xlim([0, 3 * current.cepf.cutoffK])
    max_y = np.amax(
        (current.cepf.logtau + current.cepf.logtau_THEORY_std)[current.cepf.cutoffK:3 * current.cepf.cutoffK])
    min_y = np.amin(
        (current.cepf.logtau - current.cepf.logtau_THEORY_std)[current.cepf.cutoffK:3 * current.cepf.cutoffK])
    axis.set_ylim([min_y * 0.8, max_y * 1.2])
    axis.set_xlabel(r'$P^*$')
    axis.set_ylabel(r'$L_0(P*)$')
    return axis


def plot_kappa_Pstar(current, *, axis=None, label=None, FIGSIZE=None, pstar_max=None, kappa_SI_min=None,
                     kappa_SI_max=None, pstar_tick=None, kappa_tick=None):
    """
    Plots the value of kappa as a function of P*.
    :param current: current object to plot
    :param axis: matplotlib.axes.Axes object (if None, create one)
    :param label:
    :param FIGSIZE: size of the plot

    :return: a matplotlib.axes.Axes object
    """
    if axis is None:
        figure, axis = plt.subplots(1, figsize=FIGSIZE)
    color = next(axis._get_lines.prop_cycler)['color']
    axis.fill_between(
        np.arange(current.NFREQS) + 1, (current.cepf.tau - current.cepf.tau_THEORY_std) * current.KAPPA_SCALE * 0.5,
        (current.cepf.tau + current.cepf.tau_THEORY_std) * current.KAPPA_SCALE * 0.5, alpha=0.3, color=color)
    axis.plot(np.arange(current.NFREQS) + 1, current.cepf.tau * current.KAPPA_SCALE * 0.5, 'o-', c=color, label=label)
    axis.axvline(x=current.cepf.aic_Kmin + 1, ls=':', c=color)
    axis.axvline(x=current.cepf.cutoffK + 1, ls='--', c=color)
    axis.axhline(y=current.kappa, ls='--', c=color)
    if pstar_max is None:
        pstar_max = int(round((current.cepf.cutoffK + 1) * 2.5))
    axis.set_xlim([0, pstar_max])
    if kappa_SI_max is None:
        kappa_SI_max = 1.2 * np.amax(current.KAPPA_SCALE * 0.5 *
                                     (current.cepf.tau + current.cepf.tau_THEORY_std)[current.cepf.cutoffK:pstar_max])
    if kappa_SI_min is None:
        kappa_SI_min = 0.8 * np.amin(current.KAPPA_SCALE * 0.5 *
                                     (current.cepf.tau - current.cepf.tau_THEORY_std)[current.cepf.cutoffK:pstar_max])
    axis.set_ylim([kappa_SI_min, kappa_SI_max])
    axis.set_xlabel(r'$P^*$')
    axis.set_ylabel(r'$\kappa(P^*)$ [{}]'.format(current._KAPPA_SI_UNITS))
    if pstar_tick is None:
        dx1, dx2 = _n_tick_in_range(0, pstar_max, 5)
    else:
        dx1, dx2 = (pstar_tick, pstar_tick / 2)
    if kappa_tick is None:
        dy1, dy2 = _n_tick_in_range(0, kappa_SI_max, 5)
    else:
        dy1, dy2 = (kappa_tick, kappa_tick / 2)
    axis.xaxis.set_major_locator(MultipleLocator(dx1))
    axis.xaxis.set_minor_locator(MultipleLocator(dx2))
    axis.yaxis.set_major_locator(MultipleLocator(dy1))
    axis.yaxis.set_minor_locator(MultipleLocator(dy2))
    return axis

def plot_cepstral_spectrum(current, *, freq_units='THz', freq_scale=1.0, axes=None, kappa_units=True, FIGSIZE=None, mode='log',
                           **plot_kwargs):   # yapf: disable
    """
    Plots the cepstral spectrum.

    :param current:         current object to plot
    :param freq_units:      'thz'  [THz]
                            'red'  [omega*DT/(2*pi)]
    :param freq_scale:      rescale red frequencies by this factor (e.g. 2 --> freq = [0, 0.25])
    :param axes:            matplotlib.axes.Axes object (if None, create one)
    :param kappa_units:     plot periodograms in units of kappa (default: True) - NB: log-psd not converted
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
        axes[0].plot(current.freqs_THz, current.cepf.psd * psd_scale, **plot_kwargs)
        axes[0].set_xlim([0., current.Nyquist_f_THz])
        if mode == 'log':
            axes[1].plot(current.freqs_THz, current.cepf.logpsd, **plot_kwargs)
            axes[1].set_xlim([0., current.Nyquist_f_THz])
            axes[1].set_xlabel(r'$f$ [THz]')
    elif freq_units == 'red':
        axes[0].plot(current.freqs / freq_scale, current.cepf.psd * psd_scale, **plot_kwargs)
        axes[0].set_xlim([0., 0.5 / freq_scale])
        if mode == 'log':
            axes[1].plot(current.freqs / freq_scale, current.cepf.logpsd, **plot_kwargs)
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


def plot_fstar_analysis(currents, FSTAR_THZ_LIST, original_current=None, *, axes=None, FIGSIZE=None, **plot_kwargs):
    """
    Plots kappa(P*) as a function of the f*.
    """
    if axes is None:
        figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
        return_axes = True
    else:
        return_axes = False
    axes[0].errorbar(FSTAR_THZ_LIST, [xff.kappa for xff in currents], yerr=[xff.kappa_std for xff in currents],
                     zorder=-1, **plot_kwargs)
    axes[1].errorbar(FSTAR_THZ_LIST, [xff.cepf.logtau_cutoffK for xff in currents],
                     yerr=[xff.cepf.logtau_std_cutoffK for xff in currents], zorder=-1, **plot_kwargs)
    axes[0].xaxis.set_ticks_position('top')
    axes[0].set_ylabel(r'PSD')
    axes[0].grid()
    axes[1].xaxis.set_ticks_position('bottom')
    axes[1].set_xlabel(r'$f$ [THz]')
    axes[1].set_ylabel(r'log(PSD)')
    axes[1].grid()
    if original_current is not None:
        ax2 = [axes[0].twinx(), axes[1].twinx()]
        plot_periodogram(original_current, axes=ax2, c='0.6')
        axes[0].set_ylabel(r'$\kappa$ [{}]'.format(original_current._KAPPA_SI_UNITS))
        axes[1].set_ylabel(r'$\kappa$ [{}]'.format(original_current._KAPPA_SI_UNITS))
        axes[0].set_zorder(ax2[0].get_zorder() + 1)
        axes[1].set_zorder(ax2[1].get_zorder() + 1)
        axes[0].set_frame_on(False)
        axes[1].set_frame_on(False)
    if return_axes:
        return currents, axes, figure
    else:
        return currents, axes


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

    from sportran.current import Current
    plot_kappa_units = isinstance(x, Current)
    if not axes:
        figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
        axes = plot_periodogram(x, PSD_FILTER_W=PSD_FILTER_W, freq_units=freq_units, axes=axes, mode=mode,
                                kappa_units=plot_kappa_units)   # this also updates x.PSD_FILTER_W
    xf.plot_periodogram(freq_units=freq_units, freq_scale=TSKIP, axes=axes, mode=mode, kappa_units=plot_kappa_units)
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
        plt.plot(j2pl.freqs_THz, j2pl.cepf.psd * j2pl.KAPPA_SCALE * 0.5, c=colors[1], zorder=1)
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
