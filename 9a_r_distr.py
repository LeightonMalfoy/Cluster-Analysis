"""
- plot the aftershock distance distribution for each set
    -- only for the first generation

- plot the rate/smoothed rate of earthquake
    to change between linear/smoothed, following #TODO

- Input file:
    stacked events catalog, generated by 9_stack_events.py

@author: Litong Huang, Peking University
"""

import numpy as np
import matplotlib.pyplot as plt
from src.EqCat import *
import src.seis_utils as seis_utils
import matplotlib.colors as co
import src.data_utils as data_utils
import scipy.signal


# =================================0==============================================
#                            dir, file, params
# ================================================================================
input_file = 'data/stack_events.mat'
plot_file = 'plots/k_r_distr_hist.png' #TODO:'plots/k_r_distr_rate_smoothed.png','plots/k_r_distr_rate.png'
dPar = {'Mc':2.5,
        'minMag':4.0, # minimum magnitude of mainshocks
        'maxMag':8.0, # maximum magnitude of mainshocks
        'k': 2,
        'HD_binsize':0.1,
        'lambda_binsize':0.08,
        # fitting coefficient
        'n':-1.35, 'c':1.0,
        'q':0.35, 'd':1.2,'gamma':0.6,
         }

# =================================1==============================================
#                            load catalog and select
# ================================================================================
dCluster = data_utils.loadmat(input_file)
a_MS_mag = dCluster['a_MS_mag']
a_AS_dist = dCluster['a_AS_dist']

### use different colors to discriminate different sets
index = np.arange(0,len(a_MS_mag),1)
norm=co.Normalize(vmin=0,vmax=len(a_MS_mag)-1)
cmap=plt.cm.RdYlGn
pointcolors = plt.cm.ScalarMappable(norm, cmap)
cols = pointcolors.to_rgba(index)

plt.figure()
ax = plt.subplot()

#index=[6] # for test, only plot one set
for i in index:
    #=====================================================================
    # haversine distance
    #=====================================================================
    dist = a_AS_dist[i]
    dist = dist[dist>0]
    mag = a_MS_mag[i]

    # plot the histogram
    HD_bins = np.arange(np.amin(np.log10(dist)), np.amax(np.log10(dist)), dPar['HD_binsize'])
    HD_dist, HD_bins = np.histogram(np.log10(dist), HD_bins)
    HD_bins = 10**HD_bins
    HD_dist, HD_bins = np.array(list(zip(*zip(HD_dist, HD_bins))))
    aBinSize = [bin*10**dPar['HD_binsize'] - bin for bin in HD_bins]
    HD_dist /= aBinSize
    HD_dist /= len(dist)
    ax.loglog(HD_bins[HD_dist>0], HD_dist[HD_dist>0], '-',c=cols[i], mfc = 'none', mew = .2, label='Histogram')

    # plot the rate
    HD_bins, HD_rate = seis_utils.eqRate(dist, 24)
    HD_rate /= len(dist)

    #ax.loglog(HD_bins, HD_rate, 'r-',c=cols[i], lw=1,  label='<m>=%.2f' % a_MS_mag[i])##TODO: nullify this line

    # plot the smoothed rate
    HD_bin_s, HD_rate_smoothed = seis_utils.eqlgRate(np.log10(dist),65)
    HD_bin_s, HD_rate_smoothed=HD_bin_s[HD_rate_smoothed > 0], HD_rate_smoothed[HD_rate_smoothed > 0]
    HD_rate_smoothed = scipy.signal.savgol_filter(HD_rate_smoothed,3,1)
    ##ax.loglog(HD_bin_s[HD_rate_smoothed>0],HD_rate_smoothed[HD_rate_smoothed>0],\
    ##          'b-', c=cols[i],lw=1.5,  label='<m>=%.2f'% a_MS_mag[i]) ##TODO: use this line

    """ 
    # plot Emily's fit
    L = HD_bins[np.argmax(HD_dist)]
    sel = HD_bins >= L
    aFit = np.zeros(len(HD_bins[sel]))
    aFit =dPar['c']*HD_bins[sel]**dPar['n']
    #ax.loglog(HD_bins[sel],aFit,'g-',lw=1.5,label='Felzer & Brodsky')

    # plot Hainzl et al'f fit
    sel = HD_bins>=10
    aFit = np.zeros(len(HD_bins[sel]))
    r=HD_bins[sel]
    R =1+(10/L)**(dPar['gamma']+1)
    beta = dPar['q']/dPar['d']*(R**(dPar['d']/(1+dPar['gamma'])))/((R**(dPar['q']/(1+dPar['gamma'])))+dPar['q']/dPar['d']-1)
    aFit = beta*dPar['d']*(r**dPar['gamma'])/((L**(dPar['gamma']+1))*((r/L)**(dPar['gamma']+1)+1)**(1+dPar['d']/(dPar['gamma']+1)))
    ax.loglog(r,aFit,'b-',lw=1.5,label='Decay of aftershock density')
    """ #

ax.legend(loc='upper right', fontsize=6)
ax.set_xlabel('haversine distance (km)')
ax.set_ylabel('Linear Density')
plt.savefig(plot_file)
print("save plot file: ",plot_file)
plt.show()
plt.clf()

