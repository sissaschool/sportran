import numpy as np
import matplotlib.pyplot as plt


def freq_THz_to_red(f, DT_FS):
   return f/1000.*DT_FS


class GraphManager():

    def __init__(self,j=None):
        self.freqs_THz = None
        self.psd = None
        self.fpsd = None
        self.Nyquist_f_THz = None
        if j is not None:
            self.initialize(j)

    def initialize(self,j):
        self.freqs_THz = j.freqs_THz
        self.Nyquist_f_THz = j.Nyquist_f_THz
        self.fpsd = j.fpsd
        self.psd = j.psd

    def plot_periodogram(self,  freq_units='thz', axis=None, **plot_kwargs):
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
    #    elif freq_units == 'red':
    #        axis.plot(self.freqs / freq_scale, self.fpsd, **plot_kwargs)
    #        axis.plot(self.freqs / freq_scale, self.flogpsd, **plot_kwargs)
    #        axis.set_xlim([0., 0.5 / freq_scale])
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
