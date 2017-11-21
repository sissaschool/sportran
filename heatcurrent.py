################################################################################
###   heatcurrent API
################################################################################

import numpy as np
import md
from md.mdsample import MDSample
import matplotlib.pyplot as plt


class HeatCurrent(MDSample):
   """
   HeatCurrent API for thermo-cepstral analysis.
   Defines a HeatCurrent object with useful tools to perform analysis.

   INPUT:
      - j            the heat current time series (N * N_COMPONENTS array)
      - units        the units of current ('metal', 'real')
      - DT_FS        MD time step [fs]
      - TEMPERATURE  average temperature [K]
      - VOLUME       simulation cell volume [A^3]
   """

   def __init__(self, j, units, DT_FS, TEMPERATURE, VOLUME):
      MDSample.__init__(self, traj=j, DT_FS=DT_FS)
      self.initialize_units(units, TEMPERATURE, VOLUME, DT_FS)
      self.compute_psd()
      self.initialize_cepstral_parameters()
      return


   def initialize_units(self, units, TEMPERATURE, VOLUME, DT_FS):
      """
      Initializes the units and define the kappa_scale.
      """
      self.units = units
      self.TEMPERATURE = TEMPERATURE 
      self.VOLUME = VOLUME
      self.DT_FS = DT_FS
      if (self.units == 'metal'):
         self.kappa_scale = md.units.scale_kappa_METALtoSI(TEMPERATURE, VOLUME, 1.0) # timestep is already included in the PSD definition
      elif (self.units == 'real'):
         self.kappa_scale = md.units.scale_kappa_REALtoSI(TEMPERATURE, VOLUME, 1.0)
      else:
         raise ValueError('Units not supported.')
      return


   def initialize_cepstral_parameters(self):
      """
      Defines the parameters of the theoretical distribution of the cepstrum.
      """
      self.ck_THEORY_var, self.psd_THEORY_mean = \
          md.cepstral.multicomp_cepstral_parameters(self.Nfreqs, self.N_COMPONENTS)
      return


   def plot_periodogram(self, PSD_FILTER_W=None, freq_units='thz', freq_scale=1.0, axes=None, FIGSIZE=None):
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
       self.compute_psd()
       if PSD_FILTER_W is None:
          if self.FILTER_WINDOW_WIDTH is None:
             self.filter_psd(0.)
       else:
          if (freq_units == 'thz') or (freq_units == 'THz'):
             self.filter_psd(PSD_FILTER_W/1000.*self.DT_FS)
          elif (freq_units == 'red'):
             self.filter_psd(PSD_FILTER_W)
          else:
             raise ValueError('Units not valid.')

       if axes is None:
          figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
       plt.subplots_adjust(hspace = 0.1)
       if (freq_units == 'thz') or (freq_units == 'THz'):
          axes[0].plot(self.freqs_THz, self.fpsd)
          axes[1].plot(self.freqs_THz, self.flogpsd)
          axes[0].set_xlim([0., self.Nyquist_f_THz])
          axes[1].set_xlim([0., self.Nyquist_f_THz])
       elif (freq_units == 'red'):
          axes[0].plot(self.freqs/freq_scale, self.fpsd)
          axes[1].plot(self.freqs/freq_scale, self.flogpsd)
          axes[0].set_xlim([0., 0.5/freq_scale])
          axes[1].set_xlim([0., 0.5/freq_scale])
       else:
          raise ValueError('Units not valid.')
       axes[0].xaxis.set_ticks_position('top')
       axes[0].set_ylabel('PSD')
       axes[0].grid()
       axes[1].xaxis.set_ticks_position('bottom')
       axes[1].set_xlabel('f [THz]')
       axes[1].set_ylabel('log(PSD)')
       axes[1].grid()
       return axes


   def cepstral_analysis(self, aic_type='aic', Kmin_corrfactor=1.0):
      """
      Performs Cepstral Analysis on the heat current trajectory.
         aic_type      = the Akaike Information Criterion function used to choose the cutoff ('aic', 'aicc')
         Kmin_corrfactor = correction factor multiplied by the AIC cutoff (cutoff = Kmin_corrfactor * aic_Kmin)

      Resulting conductivity:
         ( kappa_Kmin  +/-  kappa_Kmin_std ) W/(m*K)
      """
      
      self.dct = md.CosFilter(self.logpsd, ck_theory_var=self.ck_THEORY_var, \
          psd_theory_mean=self.psd_THEORY_mean, aic_type=aic_type, Kmin_corrfactor=Kmin_corrfactor)
      self.dct.scan_filter_tau()
      self.kappa_Kmin     = self.dct.tau_Kmin     * self.kappa_scale * 0.5
      self.kappa_Kmin_std = self.dct.tau_std_Kmin * self.kappa_scale * 0.5
      print '  L_0*   = {:15f} +/- {:10f}'.format(self.dct.logtau_Kmin, self.dct.logtau_std_Kmin)
      print '  S_0*   = {:15f} +/- {:10f}'.format(self.dct.tau_Kmin, self.dct.tau_std_Kmin)
      print '-------------------------------------------------'
      print '  kappa* = {:15f} +/- {:10f}  W/mK'.format(self.kappa_Kmin, self.kappa_Kmin_std)
      print '-------------------------------------------------'
      return

################################################################################


def resample_current(x, TSKIP=None, fstar_THz=None, FILTER_W=None, plot=True, PSD_FILTER_W=None, freq_units='thz', FIGSIZE=None):
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
    if not isinstance(x, HeatCurrent):
       raise ValueError('x must be a HeatCurrent object.')
    if (TSKIP is not None) and (fstar_THz is not None):
       raise ValueError('Please specify either TSKIP or fstar_THz.')
    if TSKIP is None:
       if fstar_THz is None:
          raise ValueError('Please specify either TSKIP or fstar_THz.')
       else:
          TSKIP = int(round(x.Nyquist_f_THz/fstar_THz))
    if plot:
       figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
       axes = x.plot_periodogram(PSD_FILTER_W, freq_units, 1.0, axes)
    fstar_THz = x.Nyquist_f_THz / TSKIP
    fstar_idx = np.argmin(x.freqs_THz < fstar_THz)

    # filter and sample
    if FILTER_W is None:
       FILTER_W = TSKIP
    trajf = md.tools.filter_and_sample(x.traj, FILTER_W, TSKIP, 'rectangular')
    xf = HeatCurrent(trajf, x.units, x.DT_FS*TSKIP, x.TEMPERATURE, x.VOLUME)
    if plot:
       if (freq_units == 'thz') or (freq_units == 'THz'):
          xf.plot_periodogram(x.FILTER_WINDOW_WIDTH*1000./x.DT_FS, freq_units, TSKIP, axes)
       elif (freq_units == 'red'):
          print PSD_FILTER_W
          print x.FILTER_WINDOW_WIDTH
          xf.plot_periodogram(x.FILTER_WINDOW_WIDTH*TSKIP, freq_units, TSKIP, axes)

    print 'Original Nyquist freq  f_Ny =  {:12.5f} THz'.format(x.Nyquist_f_THz)
    print 'Resampling freq          f* =  {:12.5f} THz'.format(fstar_THz)
    print 'Sampling time         TSKIP =  {:12d} steps'.format(TSKIP)
    print '                            =  {:12.3f} fs'.format(TSKIP * x.DT_FS)
    print 'Original n. of frequencies  =  {:12d}'.format(x.Nfreqs)
    print 'Resampled n. of frequencies =  {:12d}'.format(xf.Nfreqs)
    print 'PSD      @cutoff  (pre-filter) = {:12.5f}'.format(x.fpsd[fstar_idx])
    print '                 (post-filter) = {:12.5f}'.format(xf.fpsd[-1])
    print 'log(PSD) @cutoff  (pre-filter) = {:12.5f}'.format(x.flogpsd[fstar_idx])
    print '                 (post-filter) = {:12.5f}'.format(xf.flogpsd[-1])
    print 'min(PSD)          (pre-filter) = {:12.5f}'.format(x.psd_min)
    print 'min(PSD)         (post-filter) = {:12.5f}'.format(xf.psd_min)
    print '% of original PSD Power f<f* (pre-filter)  = {:5f}'.format(np.trapz(x.psd[:fstar_idx+1]) / x.psd_power * 100.)

    if plot:
       if (freq_units == 'thz') or (freq_units == 'THz'):
          axes[0].axvline(x = fstar_THz, ls='--', c='k')
          axes[1].axvline(x = fstar_THz, ls='--', c='k')
          axes[0].set_xlim([0., x.Nyquist_f_THz])
          axes[1].set_xlim([0., x.Nyquist_f_THz])
       elif (freq_units == 'red'):
          axes[0].axvline(x = 0.5/TSKIP, ls='--', c='k')
          axes[1].axvline(x = 0.5/TSKIP, ls='--', c='k')
          axes[0].set_xlim([0., 0.5/TSKIP])
          axes[1].set_xlim([0., 0.5/TSKIP])
    return xf, axes

################################################################################

