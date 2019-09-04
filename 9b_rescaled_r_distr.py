"""
1. stack the distance distribution of aftershocks whose mainshocks are in a certain range
2. fit the magnitude -- maxima (of the distance distribution) to power law :
        maxima = c*10**(sigma*mag)
3. use the relationship to rescale r to lambda
4. plot the rescaled_r (that is, lambda) distribution

-- only for the first generation
"""

import numpy as np
import matplotlib.pyplot as plt
from src.EqCat import *
import src.data_utils as data_utils
import src.seis_utils as seis_utils
import scipy.stats
import scipy.io
import scipy.signal
import matplotlib.colors as co

# =================================0==============================================
#                            dir, file, params
# ================================================================================
input_file = 'data/stack_events.mat'
plot_file = 'plots/k_lambda_rate_smoothed.png'## TODO:'plots/k_lambda_distr_hist.png','plots/k_lambda_rate.png'
dPar = {'Mc': 2.5,
        'minMag': 4.0,  # minimum magnitude of mainshocks
        'maxMag': 8.0,  # maximum magnitude of mainshocks
        'k': 2,
        'HD_binsize': 0.5,
        'lambda_binsize': 0.1,
        # rate coefficient
        'k_win':60, #
        'k_bin':60, # number of bins
        # fitting coefficient
        'n': -1.35, 'c': 1.0,
        'q': 0.35, 'd': 1.2, 'gamma': 0.6,
        }

# =================================1==============================================
#                            load catalog and select
# ================================================================================
eqcat = EqCat()

dCluster = data_utils.loadmat(input_file)
a_MS_mag = dCluster['a_MS_mag']
a_AS_dist = dCluster['a_AS_dist']
a_AS_Rm = np.zeros(len(a_MS_mag))  # maxima of distance distribution

# find the maxima of each set
for i, dist in enumerate(a_AS_dist):
    aRateBin_HD, aRate_HD = seis_utils.eqRate(dist, 2 * int(len(dist) / 30))  # ,minK=70
    aRate_HD /= len(dist)
    sel = aRateBin_HD > 0.1
    aRateBin_HD = aRateBin_HD[sel]
    aRate_HD = aRate_HD[sel]
    iRm = np.argmax(aRate_HD)
    a_AS_Rm[i] = aRateBin_HD[iRm]

# =================================3==============================================
#                            fix single power law
# ================================================================================
print(a_AS_Rm[-1])
sigma, lgC, __, __, __ = scipy.stats.linregress(a_MS_mag, np.log10(a_AS_Rm))
C = 10 ** lgC  # unit: km--> m
print("C:%.5f, sigma:%.2f" % (C, sigma))

plt.figure()
ax = plt.subplot()
index = np.arange(0, len(a_MS_mag), 1)
norm = co.Normalize(vmin=0, vmax=len(a_MS_mag) - 1)
cmap = plt.cm.RdYlGn
pointcolors = plt.cm.ScalarMappable(norm, cmap)
cols = pointcolors.to_rgba(index)

#index = [5] # for test ##TODO
for i in index:
    # =====================================================================
    # calculate lambda -- rescaled_r
    # =====================================================================
    rescaled_r = []
    for dist in a_AS_dist[i]:
        rescaled_r.append(dist / 10 ** (sigma * a_MS_mag[i] + lgC))
    rescaled_r = np.array(rescaled_r)

    rescaled_r = rescaled_r[rescaled_r > 0]
    # =============================================
    # lambda - histogram
    # =============================================
    lg_rescaled_r = np.log10(rescaled_r[rescaled_r > 0])
    expbins = np.arange(int(np.min(lg_rescaled_r)) - .5*dPar['lambda_binsize'],
                        int(np.max(lg_rescaled_r)) + .5*dPar['lambda_binsize'], dPar['lambda_binsize'])
    aHist_lambda, aBins_lambda = np.histogram(lg_rescaled_r, expbins)
    aHist_lambda, aBins_lambda = np.array(list(zip(*zip(aHist_lambda, aBins_lambda))))  # cut to same length

    # cut the empty bins
    sel = aHist_lambda > 0
    aBins_lambda = aBins_lambda[sel]
    aHist_lambda = aHist_lambda[sel]
    binsizes = [10 ** (expbin + dPar['lambda_binsize']) - 10 ** expbin for expbin in aBins_lambda]
    aHist_lambda /= np.array(binsizes)
    aHist_lambda /= len(rescaled_r)
    aBins_lambda = np.array([10 ** bin for bin in aBins_lambda])

    #ax.loglog(aBins_lambda, aHist_lambda,'ko',  mfc='none', mew=.05)##TODO:'r-',  c=cols[i],

    # =============================================
    # lambda - rate
    # =============================================
    aRate_bin, aRate = seis_utils.eqRate(rescaled_r, dPar['k_win'])
    aRate /= len(rescaled_r)
    # ax.loglog(aRate_bin, aRate, 'r-',  c=cols[i], lw=1,label='<m>=%.2f' % a_MS_mag[i]) #TODO

    # =============================================
    # lambda - smoothed rate, bin in a logarithmic scale
    # =============================================
    aBin_s, aRate_s = seis_utils.eqlgRate(np.log10(rescaled_r),dPar['k_bin'])
    aBin_s, aRate_s = aBin_s[aRate_s>0], aRate_s[aRate_s>0]
    #ax.loglog(aBin_s, aRate_s, 'r-',  c=cols[i], lw=1,label='<m>=%.2f' % a_MS_mag[i]) #TODO
    aRate_smoothed = scipy.signal.savgol_filter(aRate_s,3,1)
    ax.loglog(aBin_s, aRate_smoothed, 'r-', c=cols[i], lw=1, label='<m>=%.2f' % a_MS_mag[i])  # TODO

ax.legend(loc='upper right', fontsize=6)
ax.set_xlabel('λ')
ax.set_ylabel('P(λ)')

plt.savefig(plot_file)
print("save plot file: ",plot_file)

plt.show()
plt.clf()
