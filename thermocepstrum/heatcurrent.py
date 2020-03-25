################################################################################
###   heatcurrent API
################################################################################

import numpy as np
from . import md
from .md.mdsample import MDSample

#import matplotlib.pyplot as plt
from thermocepstrum.utils.loadAfterPlt import plt
from thermocepstrum.utils.utils import PrintMethod
log = PrintMethod()

try:
    plt
except:
    log.write_log('Warning: plt undefined')


def freq_THz_to_red(f, DT_FS):
    return f / 1000. * DT_FS


class HeatCurrent(MDSample):
    """
    HeatCurrent API for thermo-cepstral analysis.
    Defines a HeatCurrent object with useful tools to perform analysis.

    INPUT:
     - j             the heat current time series (N * N_COMPONENTS array)
       For a multi-component fluid use a (N_FLUID_COMPONENTS * N * N_COMPONENTS array)
     - units         the units of current ('metal', 'real')
     - DT_FS         MD time step [fs]
     - TEMPERATURE   average temperature [K]
     - VOLUME        simulation cell volume [A^3]
     - PSD_FILTER_W  PSD filter window [freq_units] (optional)
     - freq_units    frequency units   [THz or red] (optional)
    """

    def __init__(self, j, units, DT_FS, TEMPERATURE, VOLUME, PSD_FILTER_W=None, freq_units='THz'):

        # check if we have a multicomponent fluid
        j = np.array(j, dtype=float)
        if (len(j.shape) == 3):
            if (j.shape[0] == 1):
                self.many_currents = False
                j = np.squeeze(j, axis=0)
            else:
                self.many_currents = True
        elif (len(j.shape) <= 2):
            self.many_currents = False
        else:
            raise ValueError('Shape of j {} not valid.'.format(j.shape))

        if self.many_currents:
            log.write_log('Using multicomponent code.')
            MDSample.__init__(self, traj=j[0], DT_FS=DT_FS)
            # initialize other MDSample objects needed to make the work
            self.otherMD = []
            for js in j[1:]:
                self.otherMD.append(MDSample(traj=js, DT_FS=DT_FS))
        else:
            log.write_log('Using single component code.')
            MDSample.__init__(self, traj=j, DT_FS=DT_FS)

        self.initialize_units(units, TEMPERATURE, VOLUME, DT_FS)

        if self.traj is not None:
            if PSD_FILTER_W is None:
                if not self.many_currents:
                    self.compute_psd()
                else:
                    self.compute_kappa_multi(others=self.otherMD)
            else:
                if (freq_units == 'thz') or (freq_units == 'THz'):
                    if not self.many_currents:
                        self.compute_psd(freq_THz_to_red(PSD_FILTER_W, DT_FS))
                    else:
                        self.compute_kappa_multi(self.otherMD, freq_THz_to_red(PSD_FILTER_W, DT_FS))
                elif (freq_units == 'red'):
                    if not self.many_currents:
                        self.compute_psd(PSD_FILTER_W)
                    else:
                        self.compute_kappa_multi(self.otherMD, PSD_FILTER_W)
                else:
                    raise ValueError('Freq units not valid.')
            self.initialize_cepstral_parameters()
        else:
            log.write_log('Warning: trajectory not initialized. You should manually initialize what you need.')

        self.dct = None
        return

    def __repr__(self):
        msg = 'HeatCurrent:\n' + super(HeatCurrent, self).__repr__()
        if self.dct is not None:
            msg += self.dct.__repr__()
        return msg

    # overrides MDSample methos
    def compute_psd(self, FILTER_WINDOW_WIDTH=None, method='trajectory', DT_FS=None, average_components=True,
                    normalize=False):  # yapf: disable
        if self.many_currents:
            if self.otherMD is None:
                raise RuntimeError('self.otherMD cannot be None (wrong/missing initialization?)')
            self.compute_kappa_multi(self.otherMD, FILTER_WINDOW_WIDTH, method, DT_FS, average_components, normalize)
            return
        super(HeatCurrent, self).compute_psd(FILTER_WINDOW_WIDTH, method, DT_FS, average_components, normalize)

    @staticmethod
    def get_units_list():
        return ['metal', 'real', 'qepw', 'gpumd', 'dlpoly']

    def initialize_units(self, units, TEMPERATURE, VOLUME, DT_FS):
        """
        Initializes the units and define the kappa_scale.
        """
        self.units = units
        self.TEMPERATURE = TEMPERATURE
        self.VOLUME = VOLUME
        self.DT_FS = DT_FS
        # timestep is already included in the PSD definition, so it will be ignored here
        if (self.units == 'metal'):
            self.kappa_scale = md.units.scale_kappa_METALtoSI(TEMPERATURE, VOLUME, 1.0)
        elif (self.units == 'real'):
            self.kappa_scale = md.units.scale_kappa_REALtoSI(TEMPERATURE, VOLUME, 1.0)
        elif (self.units == 'qepw'):
            self.kappa_scale = md.units.scale_kappa_QEPWtoSI(TEMPERATURE, VOLUME, 1.0)
        elif (self.units == 'gpumd'):
            self.kappa_scale = md.units.scale_kappa_GPUMDtoSI(TEMPERATURE, VOLUME, 1.0)
        elif (self.units == 'dlpoly'):
            self.kappa_scale = md.units.scale_kappa_DLPOLYtoSI(TEMPERATURE, VOLUME, 1.0)
        else:
            raise ValueError('Units not supported.')
        return

    def initialize_cepstral_parameters(self):
        """
        Defines the parameters of the theoretical distribution of the cepstrum.
        """
        if not self.many_currents:
            self.ck_THEORY_var, self.psd_THEORY_mean = \
                md.cepstral.multicomp_cepstral_parameters(self.Nfreqs, self.N_COMPONENTS)
        else:
            if self.ndf_chi is None:
                raise RuntimeError('self.ndf_chi cannot be None.')
            self.ck_THEORY_var, self.psd_THEORY_mean = \
                md.cepstral.multicomp_cepstral_parameters(self.Nfreqs, self.ndf_chi)
        return

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
        return

    ###################################
    ###  PLOT METHODS
    ###################################

    def plot_periodogram(self, PSD_FILTER_W=None, freq_units='thz', freq_scale=1.0, axes=None, kappa_units=False,
                         FIGSIZE=None, **plot_kwargs):   # yapf: disable
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
            if not self.many_currents:
                self.compute_psd()
            else:
                if self.otherMD is None:
                    raise ValueError('self.otherMD cannot be None (missing initialization?)')
                self.compute_kappa_multi(others=self.otherMD)
        if PSD_FILTER_W is None:
            if self.FILTER_WINDOW_WIDTH is None:
                self.filter_psd(0.)
        else:
            if (freq_units == 'thz') or (freq_units == 'THz'):
                self.filter_psd(freq_THz_to_red(PSD_FILTER_W, self.DT_FS))
            elif (freq_units == 'red'):
                self.filter_psd(PSD_FILTER_W)
            else:
                raise ValueError('Units not valid.')

        if kappa_units:
            # plot psd in units of kappa - the log(psd) is not converted
            psd_scale = 0.5 * self.kappa_scale
        else:
            psd_scale = 1.0
        if axes is None:
            figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
        plt.subplots_adjust(hspace=0.1)
        if (freq_units == 'thz') or (freq_units == 'THz'):
            axes[0].plot(self.freqs_THz, psd_scale * self.fpsd, **plot_kwargs)
            axes[1].plot(self.freqs_THz, self.flogpsd, **plot_kwargs)
            axes[0].set_xlim([0., self.Nyquist_f_THz])
            axes[1].set_xlim([0., self.Nyquist_f_THz])
        elif (freq_units == 'red'):
            axes[0].plot(self.freqs / freq_scale, psd_scale * self.fpsd, **plot_kwargs)
            axes[1].plot(self.freqs / freq_scale, self.flogpsd, **plot_kwargs)
            axes[0].set_xlim([0., 0.5 / freq_scale])
            axes[1].set_xlim([0., 0.5 / freq_scale])
        else:
            raise ValueError('Units not valid.')
        axes[0].xaxis.set_ticks_position('top')
        if kappa_units:
            axes[0].set_ylabel(r'PSD [W/mK]')
        else:
            axes[0].set_ylabel(r'PSD')
        axes[0].grid()
        axes[1].xaxis.set_ticks_position('bottom')
        axes[1].set_xlabel(r'$f$ [THz]')
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

    def plot_cepstral_spectrum(self, freq_units='thz', freq_scale=1.0, axes=None, kappa_units=True, FIGSIZE=None,
                               **plot_kwargs):   # yapf: disable
        if axes is None:
            figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
        plt.subplots_adjust(hspace=0.1)
        if kappa_units:
            psd_scale = 0.5 * self.kappa_scale
        else:
            psd_scale = 1.0
        if (freq_units == 'thz') or (freq_units == 'THz'):
            axes[0].plot(self.freqs_THz, self.dct.psd * psd_scale, **plot_kwargs)
            axes[1].plot(self.freqs_THz, self.dct.logpsd, **plot_kwargs)
            axes[0].set_xlim([0., self.Nyquist_f_THz])
            axes[1].set_xlim([0., self.Nyquist_f_THz])
        elif (freq_units == 'red'):
            axes[0].plot(self.freqs / freq_scale, self.dct.psd * psd_scale, **plot_kwargs)
            axes[1].plot(self.freqs / freq_scale, self.dct.logpsd, **plot_kwargs)
            axes[0].set_xlim([0., 0.5 / freq_scale])
            axes[1].set_xlim([0., 0.5 / freq_scale])
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
        axes[1].set_xlabel(r'$f$ [THz]')
        axes[1].set_ylabel(r'log(PSD)')
        axes[1].grid()
        return axes

    def resample_current(self, TSKIP=None, fstar_THz=None, FILTER_W=None, plot=True, PSD_FILTER_W=None,
                         freq_units='thz', FIGSIZE=None):   # yapf: disable
        return resample_current(self, TSKIP=TSKIP, fstar_THz=fstar_THz, FILTER_W=FILTER_W, plot=plot,
                                PSD_FILTER_W=PSD_FILTER_W, freq_units=freq_units)   # yapf: disable


#is this function needed?
#   def compute_kappa_multi(self, others, FILTER_WINDOW_WIDTH=None):
#      """Multi-component kappa calculation."""
#      log.write_log("HeatCurrent.compute_kappa_multi")
#      if FILTER_WINDOW_WIDTH is None:
#         FILTER_WINDOW_WIDTH = self.FILTER_WINDOW_WIDTH
#      multi_mdsample = super(HeatCurrent, self).compute_kappa_multi(others, FILTER_WINDOW_WIDTH, DT_FS=self.DT_FS, call_other=True)
#      multi_hc = HeatCurrent(None, self.units, self.DT_FS, self.TEMPERATURE, self.VOLUME, FILTER_WINDOW_WIDTH, 'red')
#      ### CHECK THAT PSD IS NOT INITIALIZED by the HC constructor
#      ### We'd actually need a HC constructor from a MDSample object
#      ### NEED TO CHANGE FILTER_WINDOW_WIDTH into PSD_FILTER_W and check freq_units (for now = 'red')
#      multi_hc.initialize_psd(psd=multi_mdsample.psd, freqs=multi_mdsample.freqs)
#      multi_hc.filter_psd(FILTER_WINDOW_WIDTH)
#      multi_hc.covarALL = multi_mdsample.covarALL
#      multi_hc.ndf_chi = multi_mdsample.ndf_chi
#      multi_hc.cospectrum = multi_mdsample.cospectrum
#      return multi_hc

################################################################################


def resample_current(x, TSKIP=None, fstar_THz=None, FILTER_W=None, plot=True, PSD_FILTER_W=None, freq_units='thz',
                     FIGSIZE=None, verbose=True):   # yapf: disable
    """
    Simulate the resampling of x.

    Parameters
    ----------
    TSKIP        = sampling time [steps]
    fstar_THz    = target cutoff frequency [THz]
    FILTER_W     = pre-sampling filter window width [steps]
    plot         = plot the PSD (True/False)
    PSD_FILTER_W = PSD filtering window width [chosen frequency units]
    freq_units   = 'thz'  THz
                   'red'  omega*DT/(2*pi)
    FIGSIZE      = plot figure size

    Returns
    -------
    xf : HeatCurrent object
        a filtered & resampled HeatCurrent
    ax : array_like, optional (if plot=True)
        an array of plot axes
    """

    if not isinstance(x, HeatCurrent):
        raise ValueError('x must be a HeatCurrent object.')
    if (TSKIP is not None) and (fstar_THz is not None):
        raise ValueError('Please specify either TSKIP or fstar_THz.')
    if TSKIP is None:
        if fstar_THz is None:
            raise ValueError('Please specify either TSKIP or fstar_THz.')
        else:
            TSKIP = int(round(x.Nyquist_f_THz / fstar_THz))
    if plot:
        figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
        axes = x.plot_periodogram(PSD_FILTER_W, freq_units, 1.0, axes=axes)   # this also updates x.FILTER_WINDOW_WIDTH

    fstar_THz = x.Nyquist_f_THz / TSKIP
    fstar_idx = np.argmin(x.freqs_THz < fstar_THz)

    # filter and sample
    if FILTER_W is None:
        FILTER_W = TSKIP
    trajf = md.tools.filter_and_sample(x.traj, FILTER_W, TSKIP, 'rectangular')

    #TODO: document the units of frequency used everywere in the library and use a consistent scheme.

    # resample filtering window width in order to use the same filtering frequency window in the plot
    # if PSD_FILTER_W was specified, then x.FILTER_WINDOW_WIDTH was updated by the previous plot function
    #if x.FILTER_WINDOW_WIDTH is not None:
    #    PSD_FILTER_W = x.FILTER_WINDOW_WIDTH * TSKIP

    # define new HeatCurrent
    if not x.many_currents:
        xf = HeatCurrent(trajf, x.units, x.DT_FS * TSKIP, x.TEMPERATURE, x.VOLUME, PSD_FILTER_W, freq_units)
    else:
        if x.otherMD is None:
            raise RuntimeError('x.otherMD cannot be none (wrong/missing initialization?)')
        # filter_and_sample also other trajectories
        yf = []
        yf.append(trajf)
        for y in x.otherMD:
            tmp = md.tools.filter_and_sample(y.traj, FILTER_W, TSKIP, 'rectangular')
            yf.append(tmp)
        xf = HeatCurrent(yf, x.units, x.DT_FS * TSKIP, x.TEMPERATURE, x.VOLUME, PSD_FILTER_W, freq_units)
    if plot:
        if (freq_units == 'thz') or (freq_units == 'THz'):
            xf.plot_periodogram(x.FILTER_WINDOW_WIDTH * 1000. / x.DT_FS, 'thz', TSKIP, axes=axes)
        elif (freq_units == 'red'):
            log.write_log(PSD_FILTER_W)
            log.write_log(x.FILTER_WINDOW_WIDTH)
            xf.plot_periodogram(x.FILTER_WINDOW_WIDTH * TSKIP, 'red', TSKIP, axes=axes)

    # write log
    xf.resample_log = \
        '-----------------------------------------------------\n' +\
        '  RESAMPLE TIME SERIES\n' +\
        '-----------------------------------------------------\n' +\
        ' Original Nyquist freq  f_Ny =  {:12.5f} THz\n'.format(x.Nyquist_f_THz) +\
        ' Resampling freq          f* =  {:12.5f} THz\n'.format(fstar_THz) +\
        ' Sampling time         TSKIP =  {:12d} steps\n'.format(TSKIP) +\
        '                             =  {:12.3f} fs\n'.format(TSKIP * x.DT_FS) +\
        ' Original  n. of frequencies =  {:12d}\n'.format(x.Nfreqs) +\
        ' Resampled n. of frequencies =  {:12d}\n'.format(xf.Nfreqs)
    if x.fpsd is not None and xf.fpsd is not None:
        xf.resample_log += \
            ' PSD      @cutoff  (pre-filter) = {:12.5f}\n'.format(x.fpsd[fstar_idx]) +\
            '                  (post-filter) = {:12.5f}\n'.format(xf.fpsd[-1]) +\
            ' log(PSD) @cutoff  (pre-filter) = {:12.5f}\n'.format(x.flogpsd[fstar_idx]) +\
            '                  (post-filter) = {:12.5f}\n'.format(xf.flogpsd[-1]) +\
            ' min(PSD)          (pre-filter) = {:12.5f}\n'.format(x.psd_min) +\
            ' min(PSD)         (post-filter) = {:12.5f}\n'.format(xf.psd_min) +\
            ' % of original PSD Power f<f* (pre-filter)  = {:5f}\n'.format(np.trapz(x.psd[:fstar_idx+1]) / x.psd_power * 100.)
    else:
        xf.resample_log += ' fPSD not calculated before resampling!\n '
    xf.resample_log += '-----------------------------------------------------\n'
    if verbose:
        log.write_log(xf.resample_log)

    if plot:
        if (freq_units == 'thz') or (freq_units == 'THz'):
            axes[0].axvline(x=fstar_THz, ls='--', c='k')
            axes[1].axvline(x=fstar_THz, ls='--', c='k')
            axes[0].set_xlim([0., x.Nyquist_f_THz])
            axes[1].set_xlim([0., x.Nyquist_f_THz])
        elif (freq_units == 'red'):
            axes[0].axvline(x=0.5 / TSKIP, ls='--', c='k')
            axes[1].axvline(x=0.5 / TSKIP, ls='--', c='k')
            axes[0].set_xlim([0., 0.5 / TSKIP])
            axes[1].set_xlim([0., 0.5 / TSKIP])
        return xf, axes
    else:
        return xf


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


################################################################################
