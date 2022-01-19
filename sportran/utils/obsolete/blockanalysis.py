# -*- coding: utf-8 -*-

## !! THIS MODULE IS DEPRECATED !!

################################################################################
###
###   Code to analyze heat current trajectories in blocks using the CepstralFilter
###
################################################################################

import numpy as np
import sportran as st
from scipy.stats import shapiro
from sportran.utils import log


class MDBlocks(object):

    def __init__(self, traj, DT_FS, TSKIP, FILTER_WIDTH_T, temperature, volume, units, GUI=False):
        """  TSKIP           sampling time step
             FILTER_WIDTH_T  filter window width [original time steps]"""

        if (len(traj.shape) > 1):
            #raise NotImplementedError('Sorry, only 1-Dimensional trajectories accepted. For now.')
            self.MULTI_COMPONENT = True
            self.N_EQUIV_COMPONENTS = traj.shape[1]
        else:
            self.MULTI_COMPONENT = False
        self.DT_FS = DT_FS   # traj physical time step [fs]
        self.TSKIP = TSKIP   # sampling time step [DT_FS]
        self.FILTER_WIDTH_T = FILTER_WIDTH_T   # filter window width [DT_FS]
        self.Temperature = temperature
        self.Volume = volume
        self.units = units
        if (self.units == 'real'):
            self.tau_scale = scale_kappa_REALtoSI(self.Temperature, self.Volume, self.DT_FS)
        elif (self.units == 'metal'):
            self.tau_scale = scale_kappa_METALtoSI(self.Temperature, self.Volume, self.DT_FS)
        else:
            raise ValueError('Units not valid.')
        self.GUI = GUI
        if self.GUI:
            from ipywidgets import FloatProgress
            from IPython.display import display
            global FloatProgress, display

        # filter & sample big trajectory
        self.y_big = st.md.tools.filter_and_sample(traj, self.FILTER_WIDTH_T, self.TSKIP, 'rectangular', \
                                            even_NSTEPS=True, detrend=False, drop_first=True )
        self.NYQUIST_F = 0.5 / TSKIP   # Nyquist frequency (rescaled) [\omega*DT/(2*pi)]
        self.TOT_TIME = self.y_big.shape[0]
        log.write_log(' TOT_TIME     = {:}'.format(self.y_big.shape))
        log.write_log(' NYQUIST_F    = {:10g} = {:10g} THz'.format(self.NYQUIST_F, self.NYQUIST_F / DT_FS * 1000))

    def segment_trajectory(self, BLOCK_SIZE_T):
        """  BLOCK_SIZE_T    block width [original time steps]"""

        self.BLOCK_SIZE = BLOCK_SIZE_T / self.TSKIP   # block width (in original time steps)
        if (self.BLOCK_SIZE % 2 == 1):
            self.BLOCK_SIZE = self.BLOCK_SIZE - 1
        self.N_BLOCKS = int(np.floor(self.TOT_TIME / self.BLOCK_SIZE))
        log.write_log(' BLOCK_SIZE   = {:10d}'.format(self.BLOCK_SIZE))
        log.write_log(' N_BLOCKS     = {:10d}'.format(self.N_BLOCKS))

        # define blocks from segments of y_big
        self.block = [
            st.md.MDSample(traj=self.y_big[L * self.BLOCK_SIZE:(L + 1) * self.BLOCK_SIZE]) for L in range(self.N_BLOCKS)
        ]

    def cepstral_analysis(self, aic_type='aic', Kmin_corrfactor=1.0, bayes_p=False, density_grid=None):
        """Perform the Cepstral Analysis on all blocks."""

        if self.GUI:
            progbar = FloatProgress(min=0, max=100)
            progbar.description = '0 %'
            display(progbar)

        self.BLOCK_NFREQS = self.BLOCK_SIZE / 2 + 1
        if self.MULTI_COMPONENT:
            log.write_log(' N_EQUIV_COMPONENTS = {:10d}'.format(self.N_EQUIV_COMPONENTS))
            self.ck_THEORY_var, self.psd_THEORY_mean = st.md.cepstral.multicomp_cepstral_parameters(
                self.BLOCK_NFREQS, self.N_EQUIV_COMPONENTS)
        self.bayes_p = bayes_p

        if (self.N_BLOCKS == 1):
            raise NotImplementedError('One block.')

        for L in range(self.N_BLOCKS):
            if self.MULTI_COMPONENT:
                self.block[L].compute_psd(DT=self.TSKIP, DT_FS=self.DT_FS, average_components=True)
                self.block[L].cepf = st.md.CepstralFilter(self.block[L].logpsd, \
                    ck_theory_var=self.ck_THEORY_var, psd_theory_mean=self.psd_THEORY_mean, aic_type=aic_type, Kmin_corrfactor=Kmin_corrfactor, normalization=self.BLOCK_SIZE)
            else:
                self.block[L].compute_psd(DT=self.TSKIP, DT_FS=self.DT_FS)
                self.block[L].cepf = st.md.CepstralFilter(self.block[L].logpsd, aic_type=aic_type,
                                                    Kmin_corrfactor=Kmin_corrfactor,
                                                    normalization=self.BLOCK_SIZE)   # theory_var=None

            self.block[L].cepf.scan_filter_tau()
            if self.bayes_p:
                self.block[L].cepf.compute_p_aic(method='ba')
                if density_grid is not None:
                    self.density_grid = density_grid
                    self.block[L].cepf.compute_logtau_density(method='ba', only_stats=False, density_grid=density_grid)
                else:
                    self.block[L].cepf.compute_logtau_density(method='ba', only_stats=True)
            if self.GUI:
                progbar.value = float(L + 1) / self.N_BLOCKS * 100.
                progbar.description = '%5.2f %%' % progbar.value

        if self.GUI:
            progbar.close()

        self.freqs = self.block[0].freqs

    def cepstral_analysis_kappa(self, other, aic_type='aic', Kmin_corrfactor=1.0, bayes_p=False,
                                density_grid=None):   #need also "other", a class with the charge current!
        """Perform the Cepstral Analysis on all blocks."""

        if self.GUI:
            progbar = FloatProgress(min=0, max=100)
            progbar.description = '0 %'
            display(progbar)

        self.BLOCK_NFREQS = self.BLOCK_SIZE / 2 + 1
        if self.MULTI_COMPONENT:
            log.write_log(' N_EQUIV_COMPONENTS = {:10d}'.format(self.N_EQUIV_COMPONENTS))
            self.ck_THEORY_var, self.psd_THEORY_mean = st.md.cepstral.multicomp_cepstral_parameters(
                self.BLOCK_NFREQS, self.N_EQUIV_COMPONENTS - 1)   #different number of degrees of freedom!
        self.bayes_p = bayes_p

        if (self.N_BLOCKS == 1):
            raise NotImplementedError('One block.')

        for L in range(self.N_BLOCKS):
            if self.MULTI_COMPONENT:
                self.block[L].compute_kappa(other=other.block[L], DT=self.TSKIP, DT_FS=self.DT_FS,
                                            average_components=True)   #different method call!
                self.block[L].cepf = st.md.CepstralFilter(self.block[L].logpsd, \
                    ck_theory_var=self.ck_THEORY_var, psd_theory_mean=self.psd_THEORY_mean, aic_type=aic_type, Kmin_corrfactor=Kmin_corrfactor)#, normalization=self.BLOCK_SIZE) #removed (personal comunication with Loris)
            else:
                self.block[L].compute_kappa(other=other.block[L], DT=self.TSKIP,
                                            DT_FS=self.DT_FS)   #different method call!
                self.block[L].cepf = st.md.CepstralFilter(self.block[L].logpsd, aic_type=aic_type,
                                                    Kmin_corrfactor=Kmin_corrfactor)

            self.block[L].cepf.scan_filter_tau()
            if self.bayes_p:
                self.block[L].cepf.compute_p_aic(method='ba')
                if density_grid is not None:
                    self.density_grid = density_grid
                    self.block[L].cepf.compute_logtau_density(method='ba', only_stats=False, density_grid=density_grid)
                else:
                    self.block[L].cepf.compute_logtau_density(method='ba', only_stats=True)
            if self.GUI:
                progbar.value = float(L + 1) / self.N_BLOCKS * 100.
                progbar.description = '%5.2f %%' % progbar.value

        if self.GUI:
            progbar.close()

        self.freqs = self.block[0].freqs

    def spsd(self):
        """Sample PSD custom generator function."""
        for col in np.transpose([blk.psd for blk in self.block]):
            yield col

    def slogpsd(self):
        """Sample log(PSD) custom generator function."""
        for col in np.transpose([blk.logpsd for blk in self.block]):
            yield col

    def logpsdK(self):
        """DCT coefficients of log(PSD) custom generator function."""
        for col in np.transpose([blk.cepf.logpsdK for blk in self.block]):
            yield col

    def flogtau(self):
        """DCT-Filtered (variable K) log(TAU) custom generator function."""
        for col in np.transpose([blk.cepf.logtau for blk in self.block]):
            yield col

    def ftau(self):
        """DCT-Filtered (variable K) TAU custom generator function."""
        for col in np.transpose([blk.cepf.tau for blk in self.block]):
            yield col

    def aic(self):
        """DCT AIC custom generator function."""
        for col in np.transpose([blk.cepf.aic for blk in self.block]):
            yield col

    def aic_Kmin(self):
        """DCT AIC_Kmin custom generator function."""
        for blk in self.block:
            yield blk.cepf.aic_Kmin

    def flogpsd_Kmin(self):
        """DCT-Filtered@aic_Kmin log(PSD) custom generator function."""
        for col in np.transpose([blk.cepf.logpsd for blk in self.block]):
            yield col

    def fpsd_Kmin(self):
        """DCT-Filtered@aic_Kmin PSD custom generator function."""
        for col in np.transpose([blk.cepf.psd for blk in self.block]):
            yield col

    def flogtau_Kmin(self):
        """DCT-Filtered@aic_Kmin log(TAU) custom generator function."""
        for blk in self.block:
            yield blk.cepf.logtau_Kmin

    def ftau_Kmin(self):
        """DCT-Filtered@aic_Kmin TAU custom generator function."""
        for blk in self.block:
            yield blk.cepf.tau_Kmin

    def FTAU_Kmin(self):
        """DCT-Filtered@aic_Kmin KAPPA custom generator function."""
        for blk in self.block:
            yield blk.cepf.tau_Kmin * 0.5 * self.tau_scale

    def flogtau_Kmin_THEORYvar(self):
        """DCT-Filtered@aic_Kmin THEORYvar(log(TAU)) custom generator function."""
        for blk in self.block:
            yield blk.cepf.logtau_var_Kmin

    def ftau_Kmin_THEORYvar(self):
        """DCT-Filtered@aic_Kmin THEORYvar(TAU) custom generator function."""
        for blk in self.block:
            yield blk.cepf.tau_var_Kmin

    def flogtau_Kmin_THEORYstd(self):
        """DCT-Filtered@aic_Kmin THEORYstd(log(TAU)) custom generator function."""
        for blk in self.block:
            yield blk.cepf.logtau_std_Kmin

    def ftau_Kmin_THEORYstd(self):
        """DCT-Filtered@aic_Kmin THEORYstd(TAU) custom generator function."""
        for blk in self.block:
            yield blk.cepf.tau_std_Kmin

    def FTAU_Kmin_THEORYstd(self):
        """DCT-Filtered@aic_Kmin THEORYstd(KAPPA) custom generator function."""
        for blk in self.block:
            yield blk.cepf.tau_std_Kmin * 0.5 * self.tau_scale

    def p_aic(self):
        """DCT AIC weights distribution custom generator function."""
        for col in np.transpose([blk.cepf.p_aic for blk in self.block]):
            yield col

    def p_aic_Kave(self):
        """DCT AIC weights distribution mean custom generator function."""
        for blk in self.block:
            yield blk.cepf.p_aic_Kave

    def p_aic_Kstd(self):
        """DCT AIC weights distribution std custom generator function."""
        for blk in self.block:
            yield blk.cepf.p_aic_Kstd

    def flogtau_density(self):
        """DCT-Filtered log(TAU) Bayesian distribution custom generator function."""
        for col in np.transpose([blk.cepf.p_logtau_density for blk in self.block]):
            yield col

    def flogtau_density_xave(self):
        """DCT-Filtered log(TAU) Bayesian distribution mean custom generator function."""
        for blk in self.block:
            yield blk.cepf.p_logtau_density_xave

    def flogtau_density_xstd(self):
        """DCT-Filtered (TAU) Bayesian distribution std custom generator function."""
        for blk in self.block:
            yield blk.cepf.p_logtau_density_xstd

    def ftau_density(self):
        """DCT-Filtered TAU Bayesian distribution custom generator function."""
        for col in np.transpose([blk.cepf.p_tau_density for blk in self.block]):
            yield col

    def ftau_density_xave(self):
        """DCT-Filtered TAU Bayesian distribution mean custom generator function."""
        for blk in self.block:
            yield blk.cepf.p_tau_density_xave

    def ftau_density_xstd(self):
        """DCT-Filtered TAU Bayesian distribution std custom generator function."""
        for blk in self.block:
            yield blk.cepf.p_tau_density_xstd

    def compute_averages(self):
        """Compute all the averages."""
        ## sample psd
        self.spsd_ave = np.mean(list(self.spsd()), axis=1)
        self.spsd_std = np.std(list(self.spsd()), axis=1, ddof=1)
        log.write_log('\n   min(psd)           =  {:12g}'.format(np.min(list(self.spsd()))))

        ## sample log(psd)
        self.slogpsd_ave = np.mean(list(self.slogpsd()), axis=1)
        self.slogpsd_std = np.std(list(self.slogpsd()), axis=1, ddof=1)

        ## compute sample psd histogram
        #self.slogpsd_histogram()
        #self.spsd_histogram()

        ## log(psd), log(tau), psd, tau THEORY VARIANCES
        self.logpsdK_THEORY_var = self.block[0].cepf.logpsdK_THEORY_var.copy()
        self.logpsdK_THEORY_std = self.block[0].cepf.logpsdK_THEORY_std.copy()
        self.logtau_THEORY_var = self.block[0].cepf.logtau_THEORY_var.copy()
        self.logtau_THEORY_std = self.block[0].cepf.logtau_THEORY_std.copy()
        self.logpsd_THEORY_mean = self.block[0].cepf.logpsd_THEORY_mean.copy()

        ## DCT c_k
        self.logpsdK_ave = np.mean(list(self.logpsdK()), axis=1)
        self.logpsdK_var = np.var(list(self.logpsdK()), axis=1)
        self.logpsdK_std = np.std(list(self.logpsdK()), axis=1)

        ## DCT-Filtered log(tau), tau
        self.flogtau_ave = np.mean(list(self.flogtau()), axis=1)
        self.flogtau_std = np.std(list(self.flogtau()), axis=1, ddof=1)
        self.ftau_ave = np.mean(list(self.ftau()), axis=1)
        self.ftau_std = np.std(list(self.ftau()), axis=1, ddof=1)

        ## DCT aic_Kmin
        self.aic_Kmin_ave = np.mean(list(self.aic_Kmin()))
        self.aic_Kmin_std = np.std(list(self.aic_Kmin()), ddof=1)
        log.write_log('   max[AIC_Kmin]      =  {:12d}'.format(np.max(list(self.aic_Kmin()))))
        log.write_log('   AIC_Kmin           =  {:12.3f} +/- {:8f}'.format(self.aic_Kmin_ave, self.aic_Kmin_std))

        if self.bayes_p:
            ## DCT p_aic
            self.p_aic_KAVE = np.mean(list(self.p_aic_Kave()))
            self.p_aic_KSTD = np.std(list(self.p_aic_Kave()))
            self.avep_aic = np.mean(list(self.p_aic()), axis=1)
            self.avep_aic_KAVE, self.avep_aic_KSTD = st.md.aic.grid_statistics(np.arange(self.BLOCK_NFREQS),
                                                                               self.avep_aic)
            log.write_log('   AIC_weight_distr   =  {:12.3f} +/- {:8f}'.format(self.p_aic_KAVE, self.p_aic_KSTD))
            log.write_log('       ave_AIC_w check:  {:12.3f} +/- {:8f}'.format(self.avep_aic_KAVE, self.avep_aic_KSTD))

        ## DCT-Filtered@aic_Kmin log(psd)
        self.flogpsd_Kmin_ave = np.nanmean(list(self.flogpsd_Kmin()), axis=1)
        self.flogpsd_Kmin_std = np.nanstd(list(self.flogpsd_Kmin()), axis=1, ddof=1)

        ## DCT-Filtered@aic_Kmin psd
        self.fpsd_Kmin_ave = np.nanmean(list(self.fpsd_Kmin()), axis=1)
        self.fpsd_Kmin_std = np.nanstd(list(self.fpsd_Kmin()), axis=1, ddof=1)

        ## DCT-Filtered@aic_Kmin log(tau), tau, TAU (SI units)
        ## These below are the sample std of flogtau_Kmin.
        ## For the THEORETICAL std check self.flogtau_Kmin_THEORYstd
        self.flogtau_Kmin_ave = np.nanmean(list(self.flogtau_Kmin()))
        self.flogtau_Kmin_std = np.nanstd(list(self.flogtau_Kmin()), ddof=1)
        if self.N_BLOCKS > 2:
            self.flogtau_Kmin_SW_pvalue = shapiro(list(self.flogtau_Kmin()))[1]
        log.write_log('\n   flogtau[@AIC_Kmin] =  {:12f} +/- {:8f}  (SW p-value = {:5f})'.format(
            self.flogtau_Kmin_ave, self.flogtau_Kmin_std, self.flogtau_Kmin_SW_pvalue))

        self.ftau_Kmin_ave = np.nanmean(list(self.ftau_Kmin()))
        self.ftau_Kmin_std = self.flogtau_Kmin_std * self.ftau_Kmin_ave   # np.std(list(self.ftau_Kmin()), ddof=1)
        log.write_log('   ftau[@AIC_Kmin]    =  {:12f} +/- {:8f}'.format(self.ftau_Kmin_ave, self.ftau_Kmin_std))

        self.FTAU_Kmin_ave = self.ftau_Kmin_ave * 0.5 * self.tau_scale
        self.FTAU_Kmin_std = self.ftau_Kmin_std * 0.5 * self.tau_scale
        log.write_log('   FTAU[@AIC_Kmin]    =  {:12f} +/- {:8f}'.format(self.FTAU_Kmin_ave, self.FTAU_Kmin_std))

        ## log(tau), tau, TAU THEORY std @average AIC Kmin
        self.flogtau_THEORY_std_aveKmin = self.logtau_THEORY_std[int(round(self.aic_Kmin_ave))]
        self.ftau_THEORY_std_aveKmin = self.flogtau_THEORY_std_aveKmin * self.ftau_Kmin_ave
        self.FTAU_THEORY_std_aveKmin = self.ftau_THEORY_std_aveKmin * 0.5 * self.tau_scale
        log.write_log('\n   THEORY_STD_flogtau[@ave_AIC_Kmin] =    {:8f}  (errcheck: {:8f} +/- {:8f})'.format(
            self.flogtau_THEORY_std_aveKmin, np.nanmean(list(self.flogtau_Kmin_THEORYstd())),
            np.nanstd(list(self.flogtau_Kmin_THEORYstd()))))
        log.write_log('   THEORY_STD_ftau[@ave_AIC_Kmin]    =    {:8f}'.format(self.ftau_THEORY_std_aveKmin))
        log.write_log('   THEORY_STD_FTAU[@ave_AIC_Kmin]    =    {:8f}'.format(self.FTAU_THEORY_std_aveKmin))

        if self.bayes_p:
            ## DCT-Filtered Bayesian distribution log(tau), tau, TAU (SI units)
            self.flogtau_density_XAVE = np.mean(list(self.flogtau_density_xave()))
            self.flogtau_density_XSTD = np.std(list(self.flogtau_density_xave()))
            self.flogtau_avedensity = np.mean(list(self.flogtau_density()), axis=1)
            self.flogtau_avedensity_XAVE, self.flogtau_avedensity_XSTD = st.md.aic.grid_statistics(
                self.density_grid, self.flogtau_avedensity)
            log.write_log('\n   flogtau[AIC_w]     =  {:12f} +/- {:8f} (errcheck: {:8f} +/- {:8f})'.format(
                self.flogtau_density_XAVE, self.flogtau_density_XSTD,
                np.mean(np.mean(list(self.flogtau_density_xstd()))), np.std(list(self.flogtau_density_xstd()))))
            log.write_log('   flogtau[ave AIC_w] =  {:12f} +/- {:8f}'.format(self.flogtau_avedensity_XAVE,
                                                                             self.flogtau_avedensity_XSTD))

            self.ftau_density_XAVE = np.mean(list(self.ftau_density_xave()))
            self.ftau_density_XSTD = np.std(list(self.ftau_density_xave()))
            #self.ftau_avedensity = np.mean(list(self.ftau_density()), axis=1)
            #self.ftau_avedensity_XAVE, self.ftau_avedensity_XSTD = ta.grid_statistics(self.density_grid, self.ftau_avedensity)
            log.write_log('   ftau[AIC_w]        =  {:12f} +/- {:8f}'.format(self.ftau_density_XAVE,
                                                                             self.ftau_density_XSTD))
            #log.write_log('   ftau[ave AIC_w]    =  {:12f} +/- {:8f}'.format(self.ftau_avedensity_XAVE, self.ftau_avedensity_XSTD))

            self.FTAU_density_XAVE = self.ftau_density_XAVE * 0.5 * self.tau_scale
            self.FTAU_density_XSTD = self.ftau_density_XSTD * 0.5 * self.tau_scale
            #self.FTAU_avedensity = self.ftau_avedensity*0.5*self.tau_scale
            #self.FTAU_avedensity_XAVE = self.ftau_avedensity_XAVE*0.5*self.tau_scale
            #self.FTAU_avedensity_XSTD = self.ftau_avedensity_XSTD*0.5*self.tau_scale
            log.write_log('   FTAU[AIC_w]        =  {:12f} +/- {:8f}'.format(self.FTAU_density_XAVE,
                                                                             self.FTAU_density_XSTD))
            #log.write_log('   FTAU[ave AIC_w]    =  {:12f} +/- {:8f}'.format(self.FTAU_avedensity_XAVE, self.FTAU_avedensity_XSTD))
        log.write_log()

    def slogpsd_histogram(self, XBINW=20, NYBINS=100, YEDGE1=-8.0, YEDGE2=8.0):
        """Compute 2D histogram of sample log(PSD)."""
        self.slogpsd_h_XBINW = XBINW
        self.slogpsd_h_yedges = np.linspace(YEDGE1, YEDGE2, NYBINS + 1)
        self.slogpsd_h = np.zeros((self.BLOCK_NFREQS / XBINW, NYBINS))
        for f in range(self.BLOCK_NFREQS / XBINW):
            self.slogpsd_h[f, :], self.slogpsd_h_yedges = np.histogram(
                np.transpose(list(self.slogpsd())[f * XBINW:(f + 1) * XBINW][:]).flatten(), bins=self.slogpsd_h_yedges,
                density=True)
        self.slogpsd_h = self.slogpsd_h.T

    def spsd_histogram(self, XBINW=20, NYBINS=200, YEDGE1=0.0, YEDGE2=10.0):
        """Compute 2D histogram of sample PSD."""
        self.spsd_h_XBINW = XBINW
        self.spsd_h_yedges = np.linspace(YEDGE1, YEDGE2, NYBINS + 1)
        self.spsd_h = np.zeros((self.BLOCK_NFREQS / XBINW, NYBINS))
        for f in range(self.BLOCK_NFREQS / XBINW):
            self.spsd_h[f, :], self.spsd_h_yedges = np.histogram(
                np.transpose(list(self.spsd())[f * XBINW:(f + 1) * XBINW][:]).flatten(), bins=self.spsd_h_yedges,
                density=True)
        self.spsd_h = self.spsd_h.T
