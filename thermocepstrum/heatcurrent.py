# -*- coding: utf-8 -*-

################################################################################
###   heatcurrent API
################################################################################

import numpy as np
from . import md
from .md.mdsample import MDSample
from .md.tools.spectrum import freq_THz_to_red, freq_red_to_THz
from .md.tools.resample import filter_and_sample

from thermocepstrum.utils.loadAfterPlt import plt
from thermocepstrum.utils.utils import PrintMethod
log = PrintMethod()

try:
    plt
except:
    log.write_log('Warning: plt undefined')


class HeatCurrent(MDSample):
    """
    HeatCurrent API for thermo-cepstral analysis.
    Defines a HeatCurrent object with useful tools to perform analysis.

    INPUT:
     - traj          the heat current time series (N * N_COMPONENTS array)
       For a multi-component fluid use a (N_FLUID_COMPONENTS * N * N_COMPONENTS array)
     - units         the units of current ('metal', 'real')
     - DT_FS         MD time step [fs]
     - TEMPERATURE   average temperature [K]
     - VOLUME        simulation cell volume [A^3]
     - PSD_FILTER_W  PSD filter window [freq_units] (optional)
     - freq_units    frequency units   [THz or red] (optional)
    """

    def __init__(self, traj, UNITS, DT_FS, TEMPERATURE, VOLUME, PSD_FILTER_W=None, freq_units='THz'):
        self.initialize_currents(traj, DT_FS)
        self.initialize_units(UNITS, TEMPERATURE, VOLUME, DT_FS)
        if self.traj is not None:
            self.compute_psd(PSD_FILTER_W, freq_units)
            self.initialize_cepstral_parameters()
        else:
            log.write_log('Warning: trajectory not initialized. You should manually initialize what you need.')

        self.dct = None

    def __repr__(self):
        msg = 'HeatCurrent:\n' +\
              '  N_CURRENTS =  {}\n'.format(self.N_CURRENTS) +\
              super().__repr__()
        if self.otherMD:
            for current in self.otherMD:
                msg += current.__repr__()
        if self.dct:
            msg += self.dct.__repr__()
        return msg

    def _get_builder(self):
        """
        Get a tuple (class, builder) that can be used to build a new object with same parameters:
          TimeSeries, builder = self._get_builder()
          new_ts = TimeSeries(**builder)
        """
        if self.MANY_CURRENTS:
            traj_array = np.row_stack(([self.traj], [j.traj for j in self.otherMD]))
        else:
            traj_array = self.traj
        builder = dict(traj=traj_array, UNITS=self.UNITS, DT_FS=self.DT_FS, TEMPERATURE=self.TEMPERATURE,
                       VOLUME=self.VOLUME, PSD_FILTER_W=self.PSD_FILTER_W_THZ, freq_units='THz')
        return type(self), builder

    def initialize_currents(self, j, DT_FS):
        # check if we have a multicomponent fluid
        j = np.array(j, dtype=float)
        if (len(j.shape) == 3):
            self.N_CURRENTS = j.shape[0]
            if (self.N_CURRENTS == 1):
                self.MANY_CURRENTS = False
                j = np.squeeze(j, axis=0)
            else:
                self.MANY_CURRENTS = True
        elif (len(j.shape) <= 2):
            self.N_CURRENTS = 1
            self.MANY_CURRENTS = False
        else:
            raise ValueError('Shape of j {} not valid.'.format(j.shape))

        if self.MANY_CURRENTS:
            log.write_log('Using multicomponent code.')
            super().__init__(traj=j[0], DT_FS=DT_FS)
            # initialize other MDSample currents
            self.otherMD = [MDSample(traj=js, DT_FS=DT_FS) for js in j[1:]]
        else:
            log.write_log('Using single component code.')
            super().__init__(traj=j, DT_FS=DT_FS)
            self.otherMD = None

    def compute_psd(self, PSD_FILTER_W=None, freq_units='THz'):
        # overrides MDSample method
        """
        Compute the periodogram from the heat current time series.
        If a PSD_FILTER_W (expressed in freq_units) is known or given, the psd is also filtered.
        The PSD is multiplied by DT_FS at the end.
        """
        if self.MANY_CURRENTS:
            if self.otherMD is None:
                raise RuntimeError('self.otherMD cannot be None (wrong/missing initialization?)')
            self.compute_kappa_multi(self.otherMD, PSD_FILTER_W, freq_units)
        else:
            super().compute_psd(PSD_FILTER_W, freq_units)

    @staticmethod
    def get_units_list():
        return ['metal', 'real', 'qepw', 'gpumd', 'dlpoly']

    def initialize_units(self, UNITS, TEMPERATURE, VOLUME, DT_FS):
        """
        Initializes the units and define the kappa_scale.
        """
        self.UNITS = UNITS
        self.TEMPERATURE = TEMPERATURE
        self.VOLUME = VOLUME
        self.DT_FS = DT_FS
        # timestep is already included in the PSD definition, so it will be ignored here
        if (self.UNITS == 'metal'):
            self.kappa_scale = md.units.scale_kappa_METALtoSI(TEMPERATURE, VOLUME, 1.0)
        elif (self.UNITS == 'real'):
            self.kappa_scale = md.units.scale_kappa_REALtoSI(TEMPERATURE, VOLUME, 1.0)
        elif (self.UNITS == 'qepw'):
            self.kappa_scale = md.units.scale_kappa_QEPWtoSI(TEMPERATURE, VOLUME, 1.0)
        elif (self.UNITS == 'gpumd'):
            self.kappa_scale = md.units.scale_kappa_GPUMDtoSI(TEMPERATURE, VOLUME, 1.0)
        elif (self.UNITS == 'dlpoly'):
            self.kappa_scale = md.units.scale_kappa_DLPOLYtoSI(TEMPERATURE, VOLUME, 1.0)
        else:
            raise ValueError('Units not supported.')

    def initialize_cepstral_parameters(self):
        """
        Defines the parameters of the theoretical distribution of the cepstrum.
        """
        if not self.MANY_CURRENTS:
            self.ck_THEORY_var, self.psd_THEORY_mean = \
                md.cepstral.multicomp_cepstral_parameters(self.Nfreqs, self.N_COMPONENTS)
        else:
            if self.ndf_chi is None:
                raise RuntimeError('self.ndf_chi cannot be None.')
            self.ck_THEORY_var, self.psd_THEORY_mean = \
                md.cepstral.multicomp_cepstral_parameters(self.Nfreqs, self.ndf_chi)

    def cepstral_analysis(self, aic_type='aic', Kmin_corrfactor=1.0, K_PSD=None):
        """
        Performs Cepstral Analysis on the heat current trajectory.
           aic_type      = the Akaike Information Criterion function used to choose the cutoff ('aic', 'aicc')
           Kmin_corrfactor = correction factor multiplied by the AIC cutoff (cutoff = Kmin_corrfactor * aic_Kmin)

        Resulting conductivity:
            appa_Kmin  +/-  kappa_Kmin_std   [W/(m*K)]
        """

        self.dct = md.CosFilter(self.logpsd, ck_theory_var=self.ck_THEORY_var, \
            psd_theory_mean=self.psd_THEORY_mean, aic_type=aic_type, Kmin_corrfactor=Kmin_corrfactor)
        self.dct.scan_filter_tau(K_PSD=K_PSD)
        self.kappa_Kmin = self.dct.tau_Kmin * self.kappa_scale * 0.5
        self.kappa_Kmin_std = self.dct.tau_std_Kmin * self.kappa_scale * 0.5

        self.cepstral_log = \
              '-----------------------------------------------------\n' +\
              '  CEPSTRAL ANALYSIS\n' +\
              '-----------------------------------------------------\n' +\
              '  AIC_Kmin  = {:d}  (P* = {:d}, corr_factor = {:4f})\n'.format(self.dct.aic_Kmin, self.dct.aic_Kmin + 1, self.dct.Kmin_corrfactor) +\
              '  L_0*   = {:18f} +/- {:10f}\n'.format(self.dct.logtau_Kmin, self.dct.logtau_std_Kmin) +\
              '  S_0*   = {:18f} +/- {:10f}\n'.format(self.dct.tau_Kmin, self.dct.tau_std_Kmin) +\
              '-----------------------------------------------------\n' +\
              '  kappa* = {:18f} +/- {:10f}  W/mK\n'.format(self.kappa_Kmin, self.kappa_Kmin_std) +\
              '-----------------------------------------------------\n'
        log.write_log(self.cepstral_log)

    ###################################
    ###  PLOT METHODS
    ###################################

    def plot_periodogram(self, PSD_FILTER_W=None, freq_units='THz', freq_scale=1.0, axes=None, kappa_units=False,
                         FIGSIZE=None, **plot_kwargs):   # yapf: disable
        """
        Plot the periodogram.
          PSD_FILTER_W  = width of the filtering window
          freq_units    = 'thz'  [THz]
                          'red'  [omega*DT/(2*pi)]
          freq_scale    = rescale red frequencies by this factor (e.g. 2 --> freq = [0, 0.25])
          kappa_units   = plot periodograms in units of kappa (default: False) - NB: log-psd not converted
          axes          = matplotlib.axes.Axes object (if None, create one)
          FIGSIZE       = size of the plot

        Returns a matplotlib.axes.Axes object.
        """
        # recompute PSD if needed
        if self.psd is None:
            self.compute_psd()
        # (re)compute filtered psd, if a window has been defined
        if (PSD_FILTER_W is not None) or (self.PSD_FILTER_W is not None):
            self.filter_psd(PSD_FILTER_W, freq_units)
        else:   # use a zero-width (non-filtering) window
            self.filter_psd(0.)
        if kappa_units:   # plot psd in units of kappa - the log(psd) is not converted
            psd_scale = 0.5 * self.kappa_scale
        else:
            psd_scale = 1.0

        if axes is None:
            figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
            plt.subplots_adjust(hspace=0.1)
        if freq_units in ('THz', 'thz'):
            axes[0].plot(self.freqs_THz, psd_scale * self.fpsd, **plot_kwargs)
            axes[1].plot(self.freqs_THz, self.flogpsd, **plot_kwargs)
            axes[0].set_xlim([0., self.Nyquist_f_THz])
            axes[1].set_xlim([0., self.Nyquist_f_THz])
            axes[1].set_xlabel(r'$f$ [THz]')
        elif (freq_units == 'red'):
            axes[0].plot(self.freqs / freq_scale, psd_scale * self.fpsd, **plot_kwargs)
            axes[1].plot(self.freqs / freq_scale, self.flogpsd, **plot_kwargs)
            axes[0].set_xlim([0., 0.5 / freq_scale])
            axes[1].set_xlim([0., 0.5 / freq_scale])
            axes[1].set_xlabel(r'$f$ [$\omega$*DT/2$\pi$]')
        else:
            raise ValueError('Frequency units not valid.')
        axes[0].xaxis.set_ticks_position('top')
        if kappa_units:
            axes[0].set_ylabel(r'PSD [W/mK]')
        else:
            axes[0].set_ylabel(r'PSD')
        axes[0].grid()
        axes[1].xaxis.set_ticks_position('bottom')
        axes[1].set_ylabel(r'log(PSD)')
        axes[1].grid()
        return axes

    def plot_ck(self, axes=None, label=None, FIGSIZE=None):
        if axes is None:
            figure, axes = plt.subplots(1, figsize=FIGSIZE)
        color = next(axes._get_lines.prop_cycler)['color']
        axes.plot(self.dct.logpsdK, 'o-', c=color, label=label)
        axes.plot(self.dct.logpsdK + self.dct.logpsdK_THEORY_std, '--', c=color)
        axes.plot(self.dct.logpsdK - self.dct.logpsdK_THEORY_std, '--', c=color)
        axes.axvline(x=self.dct.aic_Kmin, ls='--', c=color)
        axes.set_xlabel(r'$k$')
        axes.set_ylabel(r'$c_k$')
        return axes

    def plot_L0_Pstar(self, axes=None, label=None, FIGSIZE=None):
        if axes is None:
            figure, axes = plt.subplots(1, figsize=FIGSIZE)
        color = next(axes._get_lines.prop_cycler)['color']
        axes.plot(np.arange(self.Nfreqs) + 1, self.dct.logtau, '.-', c=color, label=label)
        axes.plot(np.arange(self.Nfreqs) + 1, self.dct.logtau + self.dct.logtau_THEORY_std, '--', c=color)
        axes.plot(np.arange(self.Nfreqs) + 1, self.dct.logtau - self.dct.logtau_THEORY_std, '--', c=color)
        axes.axvline(x=self.dct.aic_Kmin + 1, ls='--', c=color)
        axes.set_xlim([0, 3 * self.dct.aic_Kmin])
        max_y = np.amax((self.dct.logtau + self.dct.logtau_THEORY_std)[self.dct.aic_Kmin:3 * self.dct.aic_Kmin])
        min_y = np.amin((self.dct.logtau - self.dct.logtau_THEORY_std)[self.dct.aic_Kmin:3 * self.dct.aic_Kmin])
        axes.set_ylim([min_y * 0.8, max_y * 1.2])
        axes.set_xlabel(r'$P^*$')
        axes.set_ylabel(r'$L_0(P*)$')
        return axes

    def plot_kappa_Pstar(self, axes=None, label=None, FIGSIZE=None):
        if axes is None:
            figure, axes = plt.subplots(1, figsize=FIGSIZE)
        color = next(axes._get_lines.prop_cycler)['color']
        axes.plot(np.arange(self.Nfreqs) + 1, self.dct.tau * self.kappa_scale * 0.5, '.-', c=color, label=label)
        axes.plot(np.arange(self.Nfreqs) + 1, (self.dct.tau + self.dct.tau_THEORY_std) * self.kappa_scale * 0.5,
                  '--', c=color)   # yapf: disable
        axes.plot(np.arange(self.Nfreqs) + 1, (self.dct.tau - self.dct.tau_THEORY_std) * self.kappa_scale * 0.5,
                  '--', c=color)   # yapf: disable
        axes.axvline(x=self.dct.aic_Kmin + 1, ls='--', c=color)
        axes.axhline(y=self.kappa_Kmin, ls='--', c=color)
        axes.set_xlim([0, 3 * self.dct.aic_Kmin])
        max_y = np.amax(self.kappa_scale * 0.5 *
                        (self.dct.tau + self.dct.tau_THEORY_std)[self.dct.aic_Kmin:3 * self.dct.aic_Kmin])
        min_y = np.amin(self.kappa_scale * 0.5 *
                        (self.dct.tau - self.dct.tau_THEORY_std)[self.dct.aic_Kmin:3 * self.dct.aic_Kmin])
        axes.set_ylim([min_y * 0.8, max_y * 1.2])
        axes.set_xlabel(r'$P^*$')
        axes.set_ylabel(r'$\kappa(P^*)$ [W/(m*K)]')
        return axes

    def plot_cepstral_spectrum(self, freq_units='THz', freq_scale=1.0, axes=None, kappa_units=True, FIGSIZE=None,
                               **plot_kwargs):   # yapf: disable
        if axes is None:
            figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
        plt.subplots_adjust(hspace=0.1)
        if kappa_units:
            psd_scale = 0.5 * self.kappa_scale
        else:
            psd_scale = 1.0
        if freq_units in ('THz', 'thz'):
            axes[0].plot(self.freqs_THz, self.dct.psd * psd_scale, **plot_kwargs)
            axes[1].plot(self.freqs_THz, self.dct.logpsd, **plot_kwargs)
            axes[0].set_xlim([0., self.Nyquist_f_THz])
            axes[1].set_xlim([0., self.Nyquist_f_THz])
            axes[1].set_xlabel(r'$f$ [THz]')
        elif (freq_units == 'red'):
            axes[0].plot(self.freqs / freq_scale, self.dct.psd * psd_scale, **plot_kwargs)
            axes[1].plot(self.freqs / freq_scale, self.dct.logpsd, **plot_kwargs)
            axes[0].set_xlim([0., 0.5 / freq_scale])
            axes[1].set_xlim([0., 0.5 / freq_scale])
            axes[1].set_xlabel(r'$f$ [$\omega$*DT/2$\pi$]')
        else:
            raise ValueError('Units not valid.')
        axes[0].xaxis.set_ticks_position('top')
        axes[0].set_ylabel(r'PSD')
        if kappa_units:
            axes[0].set_ylabel(r'PSD [W/mK]')
        else:
            axes[0].set_ylabel(r'PSD')
        axes[0].grid()
        axes[1].xaxis.set_ticks_position('bottom')
        axes[1].set_ylabel(r'log(PSD)')
        axes[1].grid()
        return axes

    def resample(self, TSKIP=None, fstar_THz=None, FILTER_W=None, plot=True, PSD_FILTER_W=None,
                 freq_units='THz', FIGSIZE=None, verbose=True):   # yapf: disable
        """
        Simulate the resampling of the time series.

        Parameters
        ----------
        TSKIP        = sampling time [steps]
        fstar_THz    = target cutoff frequency [THz]
        TSKIP and fstar_THZ are mutually exclusive.

        FILTER_W     = pre-sampling filter window width [steps]
        plot         = plot the PSD [True]
        PSD_FILTER_W = PSD filtering window width [chosen frequency units]
        freq_units   = 'thz'  [THz]
                       'red'  [omega*DT/(2*pi)]
        FIGSIZE      = plot figure size
        verbose      = print log [True]

        Returns
        -------
        xf : a filtered & resampled time series object
        ax : an array of plot axes, optional (if plot=True)
        """
        xf = super().resample(TSKIP, fstar_THz, FILTER_W, False, PSD_FILTER_W, freq_units, None, verbose)

        if plot:
            fstar_THz = xf.Nyquist_f_THz
            TSKIP = int(self.Nyquist_f_THz / xf.Nyquist_f_THz)

            figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
            axes = self.plot_periodogram(PSD_FILTER_W, freq_units, axes=axes)   # this also updates self.PSD_FILTER_W
            xf.plot_periodogram(freq_units=freq_units, freq_scale=TSKIP, axes=axes)
            if freq_units in ('THz', 'thz'):
                axes[0].axvline(x=fstar_THz, ls='--', c='k')
                axes[1].axvline(x=fstar_THz, ls='--', c='k')
                axes[0].set_xlim([0., self.Nyquist_f_THz])
                axes[1].set_xlim([0., self.Nyquist_f_THz])
            elif (freq_units == 'red'):
                axes[0].axvline(x=0.5 / TSKIP, ls='--', c='k')
                axes[1].axvline(x=0.5 / TSKIP, ls='--', c='k')
                axes[0].set_xlim([0., 0.5])
                axes[1].set_xlim([0., 0.5])
            return xf, axes
        return xf

    def fstar_analysis(self, TSKIP_LIST, aic_type='aic', Kmin_corrfactor=1.0, plot=True, axes=None, FIGSIZE=None,
                       verbose=False, **plot_kwargs):   # yapf: disable
        return fstar_analysis(self, TSKIP_LIST, aic_type, Kmin_corrfactor, plot, axes, FIGSIZE, verbose, **plot_kwargs)


################################################################################


def fstar_analysis(x, TSKIP_LIST, aic_type='aic', Kmin_corrfactor=1.0, plot=True, axes=None, FIGSIZE=None,
                   verbose=False, **plot_kwargs):   # yapf: disable
    """
    Perform cepstral analysis on a set of resampled time series, to study the effect of f*.
    For each TSKIP in TSKIP_LIST, the HeatCurrent x is filtered & resampled, and then cesptral-analysed.

    Parameters
    ----------
    TSKIP_LIST    = list of sampling times [steps]
    aic_type      = the Akaike Information Criterion function used to choose the cutoff ('aic', 'aicc')
    Kmin_corrfactor = correction factor multiplied by the AIC cutoff (cutoff = Kmin_corrfactor * aic_Kmin)
    plot          = plot the PSD (True/False)
    axes          = matplotlib.axes.Axes object (if None, create one)
    FIGSIZE       = plot figure size
    verbose       = verbose output (default: False)
    **plot_kwargs = other parameters passed to plot function

    Returns
    -------
    xf : array_like
        an array of HeatCurrents, corresponding the values of TSKIP_LIST
    ax : array_like, optional (if plot=True)
        array of plot axes
    fig : figure, optional (if axes=None)
        plot figure
    """

    if not isinstance(x, HeatCurrent):
        raise ValueError('x must be a HeatCurrent object.')

    xf = []
    for TSKIP in TSKIP_LIST:
        log.write_log('TSKIP =  {:d}'.format(TSKIP))
        xff = resample_current(x, TSKIP, plot=False, verbose=verbose)
        xff.cepstral_analysis(aic_type, Kmin_corrfactor)
        xf.append(xff)
    FSTAR_THZ_LIST = [xff.Nyquist_f_THz for xff in xf]

    if plot:
        if axes is None:
            figure, ax = plt.subplots(2, sharex=True, figsize=FIGSIZE)
        else:
            ax = axes
        ax[0].errorbar(FSTAR_THZ_LIST, [xff.kappa_Kmin for xff in xf], yerr=[xff.kappa_Kmin_std for xff in xf],
                       **plot_kwargs)
        ax[1].errorbar(FSTAR_THZ_LIST, [xff.dct.logtau_Kmin for xff in xf],
                       yerr=[xff.dct.logtau_std_Kmin for xff in xf], **plot_kwargs)
        # ax[0].plot(x.freqs_THz, x.fpsd,    **plot_kwargs)
        # ax[1].plot(x.freqs_THz, x.flogpsd, **plot_kwargs)
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
            x.plot_periodogram(axes=ax2, c=color)
            ax[0].set_ylabel(r'$\kappa$ [W/(m*K)]')
            ax[1].set_ylabel(r'$\kappa$ [W/(m*K)]')
            return xf, ax, figure
        else:
            return xf, ax
    else:
        return xf
