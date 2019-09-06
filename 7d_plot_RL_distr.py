'''
plot the R/L or Haversine Distance distribution of aftershock sequences of certain mainshoks
-- background events excluded

@author:Thomas Goebel, University of Memphis
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
eqCat.loadMatBin( os.path.join( 'data', file_in))

curr_Mc = dPar['aMc'][0]
klist = [30,30,30]#150
minklst = [50]#80, ,
namelist=['Northridge',]#'Baja','HectorMine','Landers'

for name,k in zip(namelist,klist):
    # load nearest neighbor distances
    NND_file = 'data/NND_%s.mat' % name
    dNND   = data_utils.loadmat(NND_file) #,  struct_as_record=True)

    #load Landers/Joshua Tree earthquake catalog
    cat_file='data/cat_%s.mat' % name
    asCat.loadMatBin(cat_file)

    # ==================================2=============================================
    #                       compute space-time-magnitude distance, histogram
    # ================================================================================
    # select only the clustering event pairs
    sel_cl = np.log10(dNND['aNND']) <= -5 #-4.7
    # HD
    HD = dNND['aHD'][sel_cl]
    print('catalog size: ', len(HD))

    # ==================================3=============================================
    #                          plot distance decay
    # ================================================================================


    fig1 = plt.figure(1)
    ax1 = plt.subplot( 111)
    ax1.set_title( name)
    for k in [20, 50, 100, 200]:
        a_r_bin, a_dens = seis_utils.eqRate( dNND['aHD'], k)
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
    # find MS magnitude
    print(eqCat.data.keys(), asCat.data.keys())
    f_MSmag = np.round( eqCat.data['Mag'][eqCat.data['N']==dNND['aEqID_p'][0]], 1)
    print( 'MS mag.', f_MSmag, round(asCat.data['Mag'][0],1))
    print( 'compare event IDs', dNND['aEqID_p'][0], asCat.data['N'][0])


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
    fig1.savefig( 'AS_distDecay_%s.png'%(name))
    plt.show()


