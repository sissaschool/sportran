import numpy as np
import matplotlib.pyplot as plt
import thermocepstrum as tc


def freq_THz_to_red(f, DT_FS):
   return f/1000.*DT_FS


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

    def plot_periodogram(self, PSD_FILTER_W=None, freq_units='thz', freq_scale=1.0, axis=None, external_object=None, **plot_kwargs):
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
        if self.psd is None:
            raise ValueError('psd is None')

        if (freq_units == 'thz') or (freq_units == 'THz'):
            axis.plot(self.freqs_THz, self.fpsd, **plot_kwargs)
    #       axis.plot(self.freqs_THz, self.flogpsd, **plot_kwargs)
            axis.set_xlim([0., self.Nyquist_f_THz])
    #       axis.set_xlim([0., self.Nyquist_f_THz])
        elif freq_units == 'red':
            axis.plot(self.freqs / freq_scale, self.fpsd, **plot_kwargs)
    #        axis.plot(self.freqs / freq_scale, self.flogpsd, **plot_kwargs)
            axis.set_xlim([0., 0.5 / freq_scale])
    #        axis.set_xlim([0., 0.5 / freq_scale])
        else:
            raise ValueError('Units not valid.')
        axis.xaxis.set_ticks_position('top')

        axis.set_ylabel(r'PSD [W/mK]')

        axis.grid()
        axis.xaxis.set_ticks_position('bottom')
        axis.set_xlabel(r'$f$ [THz]')
    #   axis.set_ylabel(r'log(PSD)')
        axis.grid()
        return axis

    @staticmethod
    def GUI_plot_periodogram(x, PSD_FILTER_W=None, freq_units='thz', freq_scale=1.0, axes=None, kappa_units=False,
                             FIGSIZE=None, **plot_kwargs):
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
        # if axes is None:
        #    figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
        # plt.subplots_adjust(hspace = 0.1)
        if (freq_units == 'thz') or (freq_units == 'THz'):
            axes.plot(x.freqs_THz, psd_scale * x.fpsd, **plot_kwargs)
            # axes[1].plot(x.freqs_THz, x.flogpsd, **plot_kwargs)
            axes.set_xlim([0., x.Nyquist_f_THz])
            # axes[1].set_xlim([0., x.Nyquist_f_THz])
        elif (freq_units == 'red'):
            axes.plot(x.freqs / freq_scale, psd_scale * x.fpsd, **plot_kwargs)
            # axes[1].plot(x.freqs/freq_scale, x.flogpsd, **plot_kwargs)
            axes.set_xlim([0., 0.5 / freq_scale])
            # axes[1].set_xlim([0., 0.5/freq_scale])
        # else:
        #    raise ValueError('Units not valid.')
        # axes[0].xaxis.set_ticks_position('top')
        # if kappa_units:
        #    axes[0].set_ylabel(r'PSD [W/mK]')
        # else:
        #    axes[0].set_ylabel(r'PSD')
        axes.grid()
        axes.xaxis.set_ticks_position('bottom')
        axes.set_xlabel(r'$f$ [THz]')
        # axes[1].set_ylabel(r'log(PSD)')
        axes.grid()
        return axes

    def resample_current(self, x, TSKIP=None, fstar_THz=None, FILTER_W=None, plot=True, PSD_FILTER_W=None, freq_units='thz',
                         FIGSIZE=None, axis=None, external_object=None):
        """
        Simulate the resampling of x.
          TSKIP        = sampling time [steps]
          fstar_THz    = target cutoff frequency [THz]
          FILTER_W     = pre-sampling filter window width [steps]
          plot         = plot the PSD (True/False)
          PSD_FILTER_W = PSD filtering window width [chosen frequency units]
          freq_units   = 'thz'  THz
                         'red'  omega*DT/(2*pi)
          FIGSIZE      = plot figure size
        """
        # if not isinstance(x, HeatCurrent):
        #     raise ValueError('x must be a HeatCurrent object.')
        # if (TSKIP is not None) and (fstar_THz is not None):
        #     raise ValueError('Please specify either TSKIP or fstar_THz.')
        if TSKIP is None:
            if fstar_THz is None:
                raise ValueError('Please specify either TSKIP or fstar_THz.')
            else:
                TSKIP = int(round(x.Nyquist_f_THz / fstar_THz))
        fstar_THz = x.Nyquist_f_THz / TSKIP
        # fstar_idx = np.argmin(x.freqs_THz < fstar_THz)

        # filter and sample
        if FILTER_W is None:
            FILTER_W = TSKIP
        trajf = tc.md.tools.filter_and_sample(x.traj, FILTER_W, TSKIP, 'rectangular')
        if not x.multicomponent:
            xf = tc.heatcurrent.HeatCurrent(trajf, x.units, x.DT_FS * TSKIP, x.TEMPERATURE, x.VOLUME, x.FILTER_WINDOW_WIDTH * TSKIP)
        else:
            if x.otherMD is None:
                raise ValueError('x.otherMD cannot be none (wrong/missing initialization?)')
            # filter_and_sample also other trajectories
            yf = []
            yf.append(trajf)
            for y in x.otherMD:
                tmp = tc.md.tools.filter_and_sample(y.traj, FILTER_W, TSKIP, 'rectangular')
                yf.append(tmp)
            xf = tc.heatcurrent.HeatCurrent(yf, x.units, x.DT_FS * TSKIP, x.TEMPERATURE, x.VOLUME, x.FILTER_WINDOW_WIDTH * TSKIP)
        if plot:
            if (freq_units == 'thz') or (freq_units == 'THz'):
                self.GUI_plot_periodogram(xf, xf.FILTER_WINDOW_WIDTH * 1000. / xf.DT_FS, 'thz', TSKIP, axes=axis)
            elif freq_units == 'red':
                self.GUI_plot_periodogram(xf, xf.FILTER_WINDOW_WIDTH * TSKIP, 'red', TSKIP, axes=axis)

        # xf.resample_log = '-----------------------------------------------------\n' + \
        #                   '  RESAMPLE TIME SERIES\n' + \
        #                   '-----------------------------------------------------\n' + \
        #                   ' Original Nyquist freq  f_Ny =  {:12.5f} THz\n'.format(x.Nyquist_f_THz) + \
        #                   ' Resampling freq          f* =  {:12.5f} THz\n'.format(fstar_THz) + \
        #                   ' Sampling time         TSKIP =  {:12d} steps\n'.format(TSKIP) + \
        #                   '                             =  {:12.3f} fs\n'.format(TSKIP * x.DT_FS) + \
        #                   ' Original  n. of frequencies =  {:12d}\n'.format(x.Nfreqs) + \
        #                   ' Resampled n. of frequencies =  {:12d}\n'.format(xf.Nfreqs) + \
        #                   ' PSD      @cutoff  (pre-filter) = {:12.5f}\n'.format(x.fpsd[fstar_idx]) + \
        #                   '                  (post-filter) = {:12.5f}\n'.format(xf.fpsd[-1]) + \
        #                   ' log(PSD) @cutoff  (pre-filter) = {:12.5f}\n'.format(x.flogpsd[fstar_idx]) + \
        #                   '                  (post-filter) = {:12.5f}\n'.format(xf.flogpsd[-1]) + \
        #                   ' min(PSD)          (pre-filter) = {:12.5f}\n'.format(x.psd_min) + \
        #                   ' min(PSD)         (post-filter) = {:12.5f}\n'.format(xf.psd_min) + \
        #                   ' % of original PSD Power f<f* (pre-filter)  = {:5f}\n'.format(
        #                       np.trapz(x.psd[:fstar_idx + 1]) / x.psd_power * 100.) + \
        #                   '-----------------------------------------------------\n'
        # print(xf.resample_log)

        if plot:
            if (freq_units == 'thz') or (freq_units == 'THz'):
                axis.axvline(x=fstar_THz, ls='--', c='k')
                axis.set_xlim([0., x.Nyquist_f_THz])
            elif (freq_units == 'red'):
                axis.axvline(x=0.5 / TSKIP, ls='--', c='k')
                axis.set_xlim([0., 0.5 / TSKIP])
        if external_object:
            external_object.xf = xf

    @staticmethod
    def plot_cepstral_spectrum(x, freq_units='thz', freq_scale=1.0, axis=None, kappa_units=True, FIGSIZE=None, external_object=None,
                               **plot_kwargs):
        if axis is None:
            figure, axis = plt.subplots(2, sharex=True, figsize=FIGSIZE)
        plt.subplots_adjust(hspace=0.1)
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
        axis.set_ylabel(r'log(PSD)')
        axis.grid()
        return axis
