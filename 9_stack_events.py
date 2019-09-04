"""
stack the distance distribution of aftershocks whose mainshocks are in a certain range
-- only for the first generation

output: mat_file:
        { 'a_MS_mag': np.array([]), # array of average mainshock magnitude of each set
          'a_AS_dist': np.array([]), #array of aftershock r distribution of each set
          }
"""

from selector.getBigEQ import getMajor
import numpy as np
from src.EqCat import *
import scipy.io
import src.data_utils as data_utils

# =================================0==============================================
#                            dir, file, params
# ================================================================================

dPar = {'Mc':2.5,
        'minMag':4.0, # minimum magnitude of mainshocks
        'maxMag':8.0, # maximum magnitude of mainshocks
        'k': 2,
        'HD_binsize':0.05,
        'lambda_binsize':0.08,
        # fitting coefficient
        'n':-1.35, 'c':1.0,
        'q':0.35, 'd':1.2,'gamma':0.6,
         }

# =================================1==============================================
#                            load catalog and select
# ================================================================================
eqcat = EqCat()

mins = np.arange(3.5,7.0,0.5) # TODO: min 2.6, 3.0, 3.5
maxs = mins + 0.5
a_MS_mag = []  # average mainshock magnitude of each set
a_AS_dist = [] # aftershock haversine distance distribution of each set

# load earthquake set, mainshocks [min,max]
for min,max in zip(mins,maxs):
    aCluster, dNND = getMajor(dPar['Mc'], min=min, max=max)
    aMag = np.zeros(len(aCluster))
    aDist = np.array([])
    sel = aMag >= 0
    for i,cluster in enumerate(aCluster):
        aMag[i] = cluster.data['Mag'][0]
        sel_c = cluster.data['sel_c']  # offspring
        aNND = dNND['aNND'][sel_c]
        aHD = dNND['aHD'][sel_c]
        sel_cl = np.log10(aNND) <= -5
        aHD = aHD[sel_cl]
        aDist = np.concatenate((aDist,aHD),axis=0)

    a_MS_mag.append(np.mean(aMag))
    a_AS_dist.append(aDist)

scipy.io.savemat('data/stack_events.mat',{'a_MS_mag':a_MS_mag,'a_AS_dist':a_AS_dist},do_compression  = True)


dCluster = data_utils.loadmat('data/stack_events.mat')
print(dCluster['a_MS_mag'])
print(dCluster['a_AS_dist'])