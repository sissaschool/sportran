# -*- coding: utf-8 -*-

import numpy as np


def dct_AIC(yk, theory_var=None):
    """AIC[K] = sum_{k>K} c_k^2/theory_var + 2*(K+1)
    Assumiamo di tenere tutti i k <= K."""
    aic = np.zeros(yk.size)
    if theory_var is None:
        N = 2 * (yk.size - 1)
        theory_var1 = np.pi**2 / 3. / N   # k = {0, N/2}
        theory_var2 = np.pi**2 / 6. / N   # otherwise
        for K in range(yk.size - 1):
            # aic[K] = np.sum(yk[K+1:]**2)/theory_var + 2.*K
            aic[K] = ((1. / theory_var2) * np.sum(yk[K + 1:-1]**2) + (1. / theory_var1) * yk[-1]**2) + 2. * (K + 1)
        aic[-1] = 2. * yk.size
    else:
        aic[-1] = 0.   # + (2*(yk.size+1))
        for K in range(yk.size - 2, -1, -1):   # N-2, N-3, ..., 0
            aic[K] = aic[K + 1] + yk[K + 1]**2 / theory_var[K + 1]
        aic = aic + 2. * (np.arange(yk.size) + 1)
    return aic


def dct_AICc(yk, theory_var=None):
    """AICc[K] = AIC[K] + 2*(K+1)*(K+2)/(NF-K-2)
    Assumiamo di tenere tutti i k <= K."""
    aic = dct_AIC(yk, theory_var)
    KK = np.arange(yk.size - 2)
    aic[:-2] = aic[:-2] + 2. * (KK + 1) * (KK + 2) / (yk.size - KK - 2.)
    aic[-2] = aic[-3]   # not defined
    aic[-1] = aic[-3]   # not defined
    return aic


def dct_aic_ab(yk, theory_var, A=1.0, B=2.0):
    """AIC[K] = sum_{k>K} c_k^2/theory_var + 2*K
    Assumiamo di tenere tutti i k <= K."""
    aic = np.zeros(yk.size)
    aic[-1] = 0.
    for K in range(yk.size - 2, -1, -1):   # N-2, N-3, ..., 0
        aic[K] = aic[K + 1] + yk[K + 1]**2 / theory_var[K + 1]
    aic = A * aic + B * (np.arange(yk.size) + 1)
    return aic


################################################################################
####  Bayesian method
################################################################################
#
# Produce_density takes as inputs the following numpy arrays:
#
#                INPUT :
#
# - aic: The vector with the estimate of the Akaike information aic(k), with k ranging from k_beg:k_max.
#         I think it is important to have values till k_max but the first k_beg can be different from one ;
# - sigma, mean : for the k_beg:k_max, in the same order, the mean and sigmas of the estimate
#                 of the transport coefficient ;
# - method : a string equal to 'one' or 'two' with two different methods of inizializing p[k], the estimate probability
#            of a given k. I hope in real cases the two methods should be equivalent.
#
#                RETURNS:
#
# The probability p[ik], a grid and a density for the probability of the transport coefficient on that grid, obtained as:
# density[igrid] \sim Sum_ik p[ik] * N[ mean[ik],sigma[ik] ].
# All pi factors are never inserted and taken care I hope by final normalization.
#
# The functions produce_p and produce_p_density distinguish the two conceptual steps of producing p and using this quantity
# to provide the final estimate on the transport coefficient
#


def produce_p(aic, method='ba', force_normalize=False):
    k0 = np.argmin(aic)
    kM = len(aic)
    p = np.zeros(kM)

    if (method == 'min'):
        # just the min(aic)
        p[k0] = 1.0

    elif (method == 'baroni'):
        for ik in range(kM):
            delta_aic = aic[ik] - aic[k0]
            GAMMA = 3.9215536345675050924567623117545   # (epsilon/sigma)^2
            p[ik] = np.exp(-0.25 * np.log(1. + GAMMA) * delta_aic)

    elif ((method == 'burnham-anderson') or (method == 'ba')):
        delta_aic = aic - aic[k0]
        p = np.exp(-0.5 * delta_aic)
        p = p / np.sum(p)

    elif ((method == 'burnham-anderson2') or (method == 'ba2')):
        delta_aic = aic - aic[k0]
        p = np.exp(-delta_aic)
        p = p / np.sum(p)

    elif (method == 'two'):
        for ik in range(kM):
            delta_aic = aic[ik] - aic[k0]
            p[ik] = np.exp(-delta_aic**2 / (2.0 * (kM - k0))) / np.sqrt(kM - k0)

    elif (method == 'four'):
        for ik in range(kM):
            delta_aic = aic[ik] - aic[k0]
            p[ik] = np.exp(-delta_aic**2 / (2.0 * np.abs(ik - k0)))
        p[k0] = 1.0

    else:
        raise KeyError('P_AIC METHOD not valid.')
    #p[ik] = np.exp(- delta_aic ** 2 / ( 2.0 * ( kM - ik) ) ) / np.sqrt(kM - ik)
    #p[ik] = np.exp(-delta_aic ** 2 / ( 2.0 * np.abs(ik - k0) )) / np.sqrt(np.abs(ik - k0))

    # normalize p
    if force_normalize:
        integral = np.trapz(p)
        for ik in range(kM):
            p[ik] = p[ik] / integral

    # checks
    if any(np.isnan(p)):
        raise Warning('WARNING: P contains NaN.')
    if (np.min(p < 0)):
        raise Warning('WARNING: P < 0.')
    return p


def produce_p_density(p, sigma, mean, grid=None, grid_size=1000):
    kM = len(mean)
    dm = np.min(mean)
    dM = np.max(mean)
    argdm = np.argmin(mean)
    argdM = np.argmax(mean)
    sigmam = sigma[argdm]
    sigmaM = sigma[argdM]
    if grid is None:
        return_grid = True
        # generate grid with grid_size points
        delta = ((dM + 5. * sigmaM) - (dm - 5. * sigmam)) / float(grid_size - 1)
        grid = np.linspace(dm - 5. * sigmam, dM + 5. * sigmaM, grid_size)
    else:
        return_grid = False
        delta = grid[1] - grid[0]
    density = np.zeros(len(grid))
    for ik in range(kM):
        density = density + p[ik] * np.exp(-(grid - mean[ik])**2 / (2.0 * (sigma[ik]**2))) / sigma[ik]
    somma = np.trapz(density) * delta
    density = density / somma
    if return_grid:
        return density, grid
    else:
        return density


def grid_statistics(grid, density, grid2=None):
    """Compute distribution mean and std.
      media   = \\sum_i (density[i] * grid[i])
      std     = sqrt( \\sum_i (density[i] * grid[i]^2) - media^2 )
       oppure = sqrt( \\sum_i (density[i] * grid2[i])  - media^2 )"""
    somma = np.sum(density)
    media = np.dot(density, grid) / somma
    if grid2 is None:
        var = np.dot(density, grid**2) / somma - media**2
    else:
        var = np.dot(density, grid2) / somma - media**2
    std = np.sqrt(var)
    return media, std
