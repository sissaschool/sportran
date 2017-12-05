################################################################################
###   heatcurrent API
################################################################################

import numpy as np
import md
from md.mdsample import MDSample
import matplotlib.pyplot as plt


def freq_THz_to_red(f, DT_FS):
   return f/1000.*DT_FS


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

   def __init__(self, j, units, DT_FS, TEMPERATURE, VOLUME, PSD_FILTER_W=None, freq_units='THz'):
      MDSample.__init__(self, traj=j, DT_FS=DT_FS)
      self.initialize_units(units, TEMPERATURE, VOLUME, DT_FS)
      if (freq_units == 'thz') or (freq_units == 'THz'):
         self.compute_psd(freq_THz_to_red(PSD_FILTER_W, DT_FS))
      elif (freq_units == 'red'):
         self.compute_psd(PSD_FILTER_W)
      else:
         raise ValueError('Units not valid.')
      self.initialize_cepstral_parameters()
      return


   def __repr__(self):
        msg = 'HeatCurrent:\n' + super(HeatCurrent, self).__repr__() \
                               + self.dct.__repr__()
        return msg


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


   def plot_periodogram(self, PSD_FILTER_W=None, freq_units='thz', freq_scale=1.0, axes=None, FIGSIZE=None, **plot_kwargs):
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
            self.filter_psd(freq_THz_to_red(PSD_FILTER_W, DT_FS))
         elif (freq_units == 'red'):
            self.filter_psd(PSD_FILTER_W)
         else:
            raise ValueError('Units not valid.')

      if axes is None:
         figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
      plt.subplots_adjust(hspace = 0.1)
      if (freq_units == 'thz') or (freq_units == 'THz'):
         axes[0].plot(self.freqs_THz, self.fpsd,    **plot_kwargs)
         axes[1].plot(self.freqs_THz, self.flogpsd, **plot_kwargs)
         axes[0].set_xlim([0., self.Nyquist_f_THz])
         axes[1].set_xlim([0., self.Nyquist_f_THz])
      elif (freq_units == 'red'):
         axes[0].plot(self.freqs/freq_scale, self.fpsd,    **plot_kwargs)
         axes[1].plot(self.freqs/freq_scale, self.flogpsd, **plot_kwargs)
         axes[0].set_xlim([0., 0.5/freq_scale])
         axes[1].set_xlim([0., 0.5/freq_scale])
      else:
         raise ValueError('Units not valid.')
      axes[0].xaxis.set_ticks_position('top')
      axes[0].set_ylabel(r'PSD')
      axes[0].grid()
      axes[1].xaxis.set_ticks_position('bottom')
      axes[1].set_xlabel(r'$f$ [THz]')
      axes[1].set_ylabel(r'log(PSD)')
      axes[1].grid()
      return axes


   def cepstral_analysis(self, aic_type='aic', Kmin_corrfactor=1.0):
      """
      Performs Cepstral Analysis on the heat current trajectory.
         aic_type      = the Akaike Information Criterion function used to choose the cutoff ('aic', 'aicc')
         Kmin_corrfactor = correction factor multiplied by the AIC cutoff (cutoff = Kmin_corrfactor * aic_Kmin)

      Resulting conductivity:
          kappa_Kmin  +/-  kappa_Kmin_std   [W/(m*K)]
      """
      
      self.dct = md.CosFilter(self.logpsd, ck_theory_var=self.ck_THEORY_var, \
          psd_theory_mean=self.psd_THEORY_mean, aic_type=aic_type, Kmin_corrfactor=Kmin_corrfactor)
      self.dct.scan_filter_tau()
      self.kappa_Kmin     = self.dct.tau_Kmin     * self.kappa_scale * 0.5
      self.kappa_Kmin_std = self.dct.tau_std_Kmin * self.kappa_scale * 0.5
      print '  AIC_Kmin  = {:d}  (P* = {:d})'.format(self.dct.aic_Kmin, self.dct.aic_Kmin + 1)
      print '  L_0*   = {:15f} +/- {:10f}'.format(self.dct.logtau_Kmin, self.dct.logtau_std_Kmin)
      print '  S_0*   = {:15f} +/- {:10f}'.format(self.dct.tau_Kmin, self.dct.tau_std_Kmin)
      print '-------------------------------------------------'
      print '  kappa* = {:15f} +/- {:10f}  W/mK'.format(self.kappa_Kmin, self.kappa_Kmin_std)
      print '-------------------------------------------------'
      return


   def plot_ck(self, axes=None, label=None, FIGSIZE=None):
      if axes is None:
         figure, axes = plt.subplots(1, figsize=FIGSIZE)
      color = next(axes._get_lines.prop_cycler)['color']
      axes.plot(self.dct.logpsdK, 'o-', c=color, label=label)
      axes.plot(self.dct.logpsdK + self.dct.logpsdK_THEORY_std, '--', c=color)
      axes.plot(self.dct.logpsdK - self.dct.logpsdK_THEORY_std, '--', c=color)
      axes.axvline(x = self.dct.aic_Kmin, ls='--', c=color)
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
      axes.axvline(x = self.dct.aic_Kmin + 1, ls='--', c=color)
      axes.set_xlim([0, 3*self.dct.aic_Kmin])
      axes.set_xlabel(r'$P^*$')
      axes.set_ylabel(r'$L_0(P*)$')
      return axes


   def plot_kappa_Pstar(self, axes=None, label=None, FIGSIZE=None):
      if axes is None:
         figure, axes = plt.subplots(1, figsize=FIGSIZE)
      color = next(axes._get_lines.prop_cycler)['color']
      axes.plot(np.arange(self.Nfreqs) + 1, self.dct.tau * self.kappa_scale * 0.5, '.-', c=color, label=label)
      axes.plot(np.arange(self.Nfreqs) + 1, (self.dct.tau + self.dct.tau_THEORY_std) * self.kappa_scale * 0.5, '--', c=color)
      axes.plot(np.arange(self.Nfreqs) + 1, (self.dct.tau - self.dct.tau_THEORY_std) * self.kappa_scale * 0.5, '--', c=color)
      axes.axvline(x = self.dct.aic_Kmin + 1, ls='--', c=color)
      axes.axhline(y = self.kappa_Kmin, ls='--', c=color)
      axes.set_xlim([0, 3*self.dct.aic_Kmin])
      axes.set_xlabel(r'$P^*$')
      axes.set_ylabel(r'$\kappa(P^*)$ [W/(m*K)]')
      return axes


   def plot_cepstral_spectrum(self, freq_units='thz', freq_scale=1.0, axes=None, FIGSIZE=None, **plot_kwargs):
       if axes is None:
          figure, axes = plt.subplots(2, sharex=True, figsize=FIGSIZE)
       plt.subplots_adjust(hspace = 0.1)
       if (freq_units == 'thz') or (freq_units == 'THz'):
          axes[0].plot(self.freqs_THz, self.dct.psd,    **plot_kwargs)
          axes[1].plot(self.freqs_THz, self.dct.logpsd, **plot_kwargs)
          axes[0].set_xlim([0., self.Nyquist_f_THz])
          axes[1].set_xlim([0., self.Nyquist_f_THz])
       elif (freq_units == 'red'):
          axes[0].plot(self.freqs/freq_scale, self.dct.psd,    **plot_kwargs)
          axes[1].plot(self.freqs/freq_scale, self.dct.logpsd, **plot_kwargs)
          axes[0].set_xlim([0., 0.5/freq_scale])
          axes[1].set_xlim([0., 0.5/freq_scale])
       else:
          raise ValueError('Units not valid.')
       axes[0].xaxis.set_ticks_position('top')
       axes[0].set_ylabel(r'PSD')
       axes[0].grid()
       axes[1].xaxis.set_ticks_position('bottom')
       axes[1].set_xlabel(r'$f$ [THz]')
       axes[1].set_ylabel(r'log(PSD)')
       axes[1].grid()
       return axes


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
   xf = HeatCurrent(trajf, x.units, x.DT_FS*TSKIP, x.TEMPERATURE, x.VOLUME, x.FILTER_WINDOW_WIDTH*TSKIP)
   if plot:
      if (freq_units == 'thz') or (freq_units == 'THz'):
         xf.plot_periodogram(x.FILTER_WINDOW_WIDTH*1000./x.DT_FS, 'thz', TSKIP, axes)
      elif (freq_units == 'red'):
         print PSD_FILTER_W
         print x.FILTER_WINDOW_WIDTH
         xf.plot_periodogram(x.FILTER_WINDOW_WIDTH*TSKIP, 'red', TSKIP, axes)

   print ' Original Nyquist freq  f_Ny =  {:12.5f} THz'.format(x.Nyquist_f_THz)
   print ' Resampling freq          f* =  {:12.5f} THz'.format(fstar_THz)
   print ' Sampling time         TSKIP =  {:12d} steps'.format(TSKIP)
   print '                             =  {:12.3f} fs'.format(TSKIP * x.DT_FS)
   print ' Original  n. of frequencies =  {:12d}'.format(x.Nfreqs)
   print ' Resampled n. of frequencies =  {:12d}'.format(xf.Nfreqs)
   print ' PSD      @cutoff  (pre-filter) = {:12.5f}'.format(x.fpsd[fstar_idx])
   print '                  (post-filter) = {:12.5f}'.format(xf.fpsd[-1])
   print ' log(PSD) @cutoff  (pre-filter) = {:12.5f}'.format(x.flogpsd[fstar_idx])
   print '                  (post-filter) = {:12.5f}'.format(xf.flogpsd[-1])
   print ' min(PSD)          (pre-filter) = {:12.5f}'.format(x.psd_min)
   print ' min(PSD)         (post-filter) = {:12.5f}'.format(xf.psd_min)
   print ' % of original PSD Power f<f* (pre-filter)  = {:5f}\n'.format(np.trapz(x.psd[:fstar_idx+1]) / x.psd_power * 100.)

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
   else:
      return xf


def fstar_analysis(x, TSKIP_LIST, aic_type='aic', Kmin_corrfactor=1.0, plot=True, axes=None, FIGSIZE=None, **plot_kwargs):
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

   xf = []
   for TSKIP in TSKIP_LIST:
      print 'TSKIP =  {:d}'.format(TSKIP)
      xff = resample_current(x, TSKIP, plot=False)
      xff.cepstral_analysis(aic_type, Kmin_corrfactor)
      xf.append( xff )
   FSTAR_THZ_LIST = [ xff.Nyquist_f_THz for xff in xf ]

   if plot:
      if axes is None:
         figure, ax = plt.subplots(2, sharex=True, figsize=FIGSIZE)
      else:
         ax = axes
      ax[0].errorbar( FSTAR_THZ_LIST, [xff.kappa_Kmin for xff in xf], yerr = [xff.kappa_Kmin_std for xff in xf], **plot_kwargs )
      ax[1].errorbar( FSTAR_THZ_LIST, [xff.dct.logtau_Kmin for xff in xf], yerr = [xff.dct.logtau_std_Kmin for xff in xf], **plot_kwargs )
      #ax[0].plot(x.freqs_THz, x.fpsd,    **plot_kwargs)
      #ax[1].plot(x.freqs_THz, x.flogpsd, **plot_kwargs)
      ax[0].xaxis.set_ticks_position('top')
      ax[0].set_ylabel(r'PSD')
      ax[0].grid()
      ax[1].xaxis.set_ticks_position('bottom')
      ax[1].set_xlabel(r'$f$ [THz]')
      ax[1].set_ylabel(r'log(PSD)')
      ax[1].grid()

      if axes is None:
         ax2 = [ax[0].twinx(), ax[1].twinx()]
         color=next(ax[0]._get_lines.prop_cycler)['color']
         color=next(ax[1]._get_lines.prop_cycler)['color']
         x.plot_periodogram(axes=ax2, c=color)
         ax[0].set_ylabel(r'$\kappa$ [W/(m*K)]')
         ax[1].set_ylabel(r'$\kappa$ [W/(m*K)]')
         return xf, ax, figure
      else:
         return xf, ax
   else:
      return xf


################################################################################

