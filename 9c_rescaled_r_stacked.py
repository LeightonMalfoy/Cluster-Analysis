'''
plot the R/L or Haversine Distance distribution of aftershock sequences of certain mainshoks
-- background events excluded
'''
import matplotlib as mpl
#mpl.use( 'Agg') # uncomment for interactive plotting
import matplotlib.pyplot as plt

import os
import numpy as np
#------------------------------my modules--------------------------------------
import src.data_utils as data_utils
import src.seis_utils as seis_utils
from src.EqCat import *
import src.plot_hist as myplot

eqCat   = EqCat() # original catalog
asCat   = EqCat()
#=================================0==============================================
#                            dir, file, params
#================================================================================
dir_in = '%s/data'%( os.path.expanduser('.'))
file_in= 'hs_1981_2018_all.mat'

#file_b  = '%s_b_Mc_D.txt'%(fileIn.split('.')[0])
dPar  = {   'aMc'         :  np.array( [2.5,2.6,3.0,3.5,4.0]),#np.array([3.0, 4.0]),np.array( [2.0, 2.5, 3.0, 3.5]
            # fractal dimension and b for eq. (1)
            'D'           : 1.6, #1.6  TODO: - these values should be constrained independently
            'b'           : 1.0,
            #=================plotting==============
            'eta_binsize' :  .2,
            'nnd_binsize' : .2,
            'hd_binsize' : .2,
            'xmin' : 0, 'xmax' : 1.0,
            'eta_0': -5.0,
          }

#=================================1==============================================
#                            load data, select events
#================================================================================
asSets = data_utils.loadmat('data/stack_events_2.5_7_atep1.mat')
a_MS_mag = asSets['a_MS_mag']
a_AS_dist = asSets['a_AS_dist']

i=0
for f_MSmag,dist in zip(a_MS_mag,a_AS_dist):

    # ==================================3=============================================
    #                          plot distance decay
    # ================================================================================

    name='<m>=%.2f'%f_MSmag
    fig1 = plt.figure(1)
    ax1 = plt.subplot( 111)
    ax1.set_title( name)
    for k in [20, 50, 100, 200]:
        a_r_bin, a_dens = seis_utils.eqRate( dist, k)
        ax1.loglog( a_r_bin, a_dens, 'o', label = str(k), mew = 0)
    #--------- plot -1.4 slope from felzer and Brodsky for comparison-----------------
    gamma = -1.4
    selPeak = a_dens == a_dens.max()
    preFac= a_dens[selPeak]/(a_r_bin[selPeak]**gamma)
    ax1.plot( a_r_bin, preFac*a_r_bin**gamma, 'w-', lw =2.5)
    ax1.plot( a_r_bin, preFac*a_r_bin**gamma, 'k--', label = 'Dens ~ r ^ %.1f'%( gamma))
    # ==================================4=============================================
    #                       highlight mainshock rupture dimension
    # ================================================================================
    l0,sigma = 0.01, 0.44
    f_L1 = l0*10**(sigma*f_MSmag) # Hainzl , Moradpour et al.
    # source dimension assuming constant stress drop of 3 MPa
    f_M0 = 10**( 1.5*( f_MSmag + 6)) # M0 in Nm

    f_deltaTau = 3*1e6# stress drop
    f_L2 = ( 4*f_M0/(np.pi*f_deltaTau))**(1./3)*1e-3
    #f_L2 = 1.99*10 ** (0.49 * asCat.data['Mag'][0]-3)

    ax1.loglog( [f_L1,f_L1], ax1.get_ylim(), 'b--', label = 'Moradpour: l/2=%.1f km'%( f_L1))
    ax1.loglog( [f_L2,f_L2], ax1.get_ylim(), 'r--', label = 'const. dtau: l/2=%.1f km'%( f_L2))
    '''
    #--------- plot Moradpour et al's fit for 10 km or more for comparison-----------------
    gamma = -1.4
    selPeak = a_dens == a_dens.max()
    preFac= a_dens[selPeak]/(a_r_bin[selPeak]**gamma)
    ax1.plot( a_r_bin, preFac*a_r_bin**gamma, 'w-', lw =2.5)
    ax1.plot( a_r_bin, preFac*a_r_bin**gamma, 'k--', label = 'Dens ~ r ^ %.1f'%( gamma))
    '''
    #-----------------------------------------------------
    ax1.legend( loc = 'lower left')
    ax1.set_xlabel( 'Distance from MS [km]')
    ax1.set_ylabel( 'Seismicity Density')
    fig1.savefig( 'plots/AS_distDecay_%d.png'%(i))
    plt.show()
    i += 1
    plt.clf()


