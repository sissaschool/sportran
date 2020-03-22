import numpy as np
import matplotlib.pyplot as plt
import thermocepstrum as tc


def freq_THz_to_red(f, DT_FS):
    return f / 1000. * DT_FS


class GraphManager:

    def __init__(self, j=None):
        self.freqs_THz = None
        self.psd = None
        self.fpsd = None
        self.freqs = None
        self.Nyquist_f_THz = None
        if j is not None:
            self.initialize(j)

    def initialize(self, j):
        self.freqs_THz = j.freqs_THz
        self.Nyquist_f_THz = j.Nyquist_f_THz
        self.fpsd = j.fpsd
        self.freqs = j.freqs
        self.psd = j.psd

    @staticmethod
    def GUI_plot_periodogram(x, PSD_FILTER_W=None, freq_units='thz', freq_scale=1.0, axis=None, kappa_units=True,
                             data=None, FIGSIZE=None, **plot_kwargs):
        """
        Plot the periodogram.
          PSD_FILTER_W  = width of the filtering window
          freq_units    = 'thz'  THz
                          'red'  omega*DT/(2*pi)
          freq_scale    = rescale red frequencies by this factor (e.g. 2 --> freq = [0, 0.25])
          axes          = matplotlib.axes.Axes object (if None, create one)
          FIGSIZE       = size of the plot

        Returns a matplotlib.axes.Axes object.
        """
        # recompute PSD if needed
        if x.psd is None:
            if not x.multicomponent:
                x.compute_psd()
            else:
                if x.otherMD is None:
                    raise ValueError('x.otherMD cannot be None (missing initialization?)')
                x.compute_kappa_multi(others=x.otherMD)
        if PSD_FILTER_W is None:
            if x.FILTER_WINDOW_WIDTH is None:
                x.filter_psd(0.)
        else:
            if (freq_units == 'thz') or (freq_units == 'THz'):
                x.filter_psd(freq_THz_to_red(PSD_FILTER_W, x.DT_FS))
            elif (freq_units == 'red'):
                x.filter_psd(PSD_FILTER_W)
            else:
                raise ValueError('Units not valid.')

        if kappa_units:
            # plot psd in units of kappa - the log(psd) is not converted
            psd_scale = 0.5 * x.kappa_scale
        else:
            psd_scale = 1.0
        # if axis is None:
        #    figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
        # plt.subplots_adjust(hspace = 0.1)
        if (freq_units == 'thz') or (freq_units == 'THz'):
            axis.plot(x.freqs_THz, psd_scale * x.fpsd, **plot_kwargs)
            # axes[1].plot(x.freqs_THz, x.flogpsd, **plot_kwargs)
            axis.set_xlim([0., x.Nyquist_f_THz])
            # axes[1].set_xlim([0., x.Nyquist_f_THz])
        elif (freq_units == 'red'):
            axis.plot(x.freqs / freq_scale, psd_scale * x.fpsd, **plot_kwargs)
            # axes[1].plot(x.freqs/freq_scale, x.flogpsd, **plot_kwargs)
            axis.set_xlim([0., 0.5 / freq_scale])
            # axes[1].set_xlim([0., 0.5/freq_scale])
        # else:
        #    raise ValueError('Units not valid.')
        # axes[0].xaxis.set_ticks_position('top')
        # if kappa_units:
        #    axes[0].set_ylabel(r'PSD [W/mK]')
        # else:
        #    axes[0].set_ylabel(r'PSD')
        axis.grid()
        axis.xaxis.set_ticks_position('bottom')
        axis.set_xlabel(r'$f$ [THz]')
        # axes[1].set_ylabel(r'log(PSD)')
        axis.grid()
        return axis

    def resample_current(self, x, TSKIP=None, fstar_THz=None, FILTER_W=None, plot=True, PSD_FILTER_W=None,
                         freq_units='thz', FIGSIZE=None, axis=None, data=None):

        if data.changes:
            print('resampling current')
            xf = x.resample_current(TSKIP=TSKIP, fstar_THz=fstar_THz, FILTER_W=FILTER_W, plot=False,
                                    PSD_FILTER_W=PSD_FILTER_W, freq_units=freq_units)
        else:
            print('not resampling current (data is same) {}'.format(data.changes))
            xf = data.xf

        if plot:
            if (freq_units == 'thz') or (freq_units == 'THz'):
                self.GUI_plot_periodogram(xf, xf.FILTER_WINDOW_WIDTH * 1000. / xf.DT_FS, 'thz', TSKIP, axis=axis)
            elif freq_units == 'red':
                self.GUI_plot_periodogram(xf, xf.FILTER_WINDOW_WIDTH * TSKIP, 'red', TSKIP, axis=axis)

        if plot:
            if (freq_units == 'thz') or (freq_units == 'THz'):
                axis.axvline(x=xf.Nyquist_f_THz, ls='--', c='k')
                axis.set_xlim([0., x.Nyquist_f_THz])
            elif (freq_units == 'red'):
                axis.axvline(x=0.5 / TSKIP, ls='--', c='k')
                axis.set_xlim([0., 0.5 / TSKIP])
        if data:
            data.xf = xf

    @staticmethod
    def plot_cepstral_spectrum(x, freq_units='thz', freq_scale=1.0, axis=None, kappa_units=True, FIGSIZE=None,
                               data=None, **plot_kwargs):
        if kappa_units:
            psd_scale = 0.5 * x.kappa_scale
        else:
            psd_scale = 1.0
        if (freq_units == 'thz') or (freq_units == 'THz'):
            axis.plot(x.freqs_THz, x.dct.psd * psd_scale, **plot_kwargs)
            # axes[1].plot(x.freqs_THz, x.dct.logpsd, **plot_kwargs)
            axis.set_xlim([0., x.Nyquist_f_THz])
            # axes[1].set_xlim([0., x.Nyquist_f_THz])
        elif (freq_units == 'red'):
            axis.plot(x.freqs / freq_scale, x.dct.psd * psd_scale, **plot_kwargs)
            # axes[1].plot(x.freqs / freq_scale, x.dct.logpsd, **plot_kwargs)
            axis.set_xlim([0., 0.5 / freq_scale])
            # axes[1].set_xlim([0., 0.5 / freq_scale])
        else:
            raise ValueError('Units not valid.')
        axis.xaxis.set_ticks_position('top')
        axis.set_ylabel(r'PSD')
        if kappa_units:
            axis.set_ylabel(r'PSD [W/mK]')
        else:
            axis.set_ylabel(r'PSD')
        axis.xaxis.set_ticks_position('bottom')
        axis.set_xlabel(r'$f$ [THz]')
        axis.grid()
        return axis
