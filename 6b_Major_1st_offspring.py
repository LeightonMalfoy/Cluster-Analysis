"""
get the calalog containing one mainshock and its offsprings (1st generation aftershocks)

output: 1. eaCat
        2. NND {'aEqID_p': np.array([]),
                'aEqID_c': np.array([]),
                'aNND': np.array([]), # Nearest Neighbor Distance
                'aHD': np.array([]), # Haversine Distance between parent and offspring

@ Litong Huang, Peking Univ
"""

from selector.getBigEQ import getMajor
import scipy.io
from src.EqCat import *


aCluster,dNND=getMajor(2.5)
name='Northridge'
cluster=aCluster[8]
cluster.saveMatBin('../data/cat_%s.mat'%name)
sel_c = cluster.data['sel_c'] # offspring
aHD = dNND['aHD'][sel_c]
NND_file = '../data/NND_%s.mat' % name
nNND = {}
nNND['aEqID_p'] = dNND['aEqID_p'][sel_c]
nNND['aEqID_c'] = dNND['aEqID_c'][sel_c]
nNND['aNND'] = dNND['aNND'][sel_c]
nNND['aHD'] = aHD
scipy.io.savemat(NND_file, nNND, do_compression=True)
