'''
Created on August 16, 2016

- compute inter-event times and distance and normalize by magnitude

@author: tgoebel
'''
import matplotlib as mpl
mpl.use( 'Agg') # uncomment for interactive plotting

import os
import numpy as np

#------------------------------my modules--------------------------------------
import src.data_utils as data_utils
import src.seis_utils as seis_utils
from src.EqCat import *
import src.plot_hist as myplot

eqCat   = EqCat() # original catalog
eqCatMc = EqCat() # this catalog wil be modfied with each Mc iteration
catChild=  EqCat()
catParent= EqCat()
#=================================1==============================================
#                            dir, file, params
#================================================================================
dir_in = '%s/data'%( os.path.expanduser('.'))
file_in= 'hs_1981_2018_all.mat'

#file_b  = '%s_b_Mc_D.txt'%(fileIn.split('.')[0])
dPar  = {   'aMc'         :  np.array( [2.5,3.0,3.5,4.0]),#np.array([3.0, 4.0]),np.array( [2.0, 2.5, 3.0, 3.5]
            # fractal dimension and b for eq. (1)
            'D'           : 1.6, #1.6  TODO: - these values should be contrained independently
            'b'           : 1.0,
            #=================plotting==============
            'eta_binsize' :  .1                                                                                                                                                      ,
            'xmin' : 0, 'xmax' : 1.0,
            'eta_0': -5.0,
            #'k' :[int(300*(i+1)**(np.log10(1/6))) for i in list(range(10))],
            'k':[600,400,200,50], # to calculate the r/l rate
          }

#=================================2==============================================
#                            load data, select events
#================================================================================
eqCat.loadMatBin(  os.path.join( dir_in, file_in))
print( 'total no. of events', eqCat.size())
eqCat.selectEvents(dPar['aMc'][0], None, 'Mag')
#eqCat.selector( tmin, tmax, 'Time')
print( 'no. of events after initial selection', eqCat.size())
# cut below current completeness
'''
curr_Mc=dPar['aMc'][2]
k=dPar['k'][2]
if True:
    eqCatMc.copy(eqCat)
    eqCatMc.selector(curr_Mc, None, 'Mag')
    print( 'current catalog size: ',eqCatMc.size())

    # load nearest neighbor distances
    NND_file = 'data/%s_NND_Mc_%.1f_HD.mat'%( file_in.split('.')[0], curr_Mc)
    dNND   = data_utils.loadmat(NND_file) #,  struct_as_record=True)
    # ==================================2=============================================
    #                       compute space-time-magnitude distance, histogram
    # ================================================================================
    ## l
    dsigma = 5e6
    Mw = eqCatMc.data['Mag'][1:]
    M0 = 10 ** (1.5 * (Mw + 6.03))
    l = 0.001*(7 / 16 * M0 / dsigma) ** (1 / 3)
    ## hd
    HD = dNND['aHD']
    ## Rl
    aRl = HD / l
    ###histogram
    aBins = np.arange(0,40, dPar['eta_binsize'])
    aHist, aBins = np.histogram(aRl, aBins)
    aHist, aBins = np.array(list(zip(*zip(aHist, aBins))))  # cut to same length
    aRate_bin, aRate = seis_utils.eqRate(aRl,k)
    aHist /= dPar['eta_binsize']
    # =================================4==============================================
    #                          plot histogram
    # ================================================================================
    #myplot.plot_rate_pd(aBins, aHist, aRate_bin, aRate, curr_Mc)
    myplot.plot_loglog(aBins,aHist,aRate_bin,aRate,curr_Mc,file_in,dPar['eta_binsize'])
'''
for curr_Mc,k in zip(dPar['aMc'],dPar['k']):
    # cut below current completeness
    eqCatMc.copy(eqCat)
    eqCatMc.selectEvents(curr_Mc, None, 'Mag')
    print( 'current catalog size: ',eqCatMc.size())

    # load nearest neighbor distances
    NND_file = 'data/%s_NND_Mc_%.1f_HD.mat'%( file_in.split('.')[0], curr_Mc)
    dNND   = data_utils.loadmat(NND_file) #,  struct_as_record=True)
    # ==================================2=============================================
    #                       compute space-time-magnitude distance, histogram
    # ================================================================================
    catChild.copy(eqCatMc)
    catParent.copy(eqCatMc)
    # catChild, catPar = create_parent_child_cat( projCat, dNND)
    catChild.selEventsFromID(dNND['aEqID_c'], repeats=True)
    catParent.selEventsFromID(dNND['aEqID_p'], repeats=True)
    print('size of parent catalog', catParent.size(), 'size of offspring cat', catChild.size())
    ## l
    dsigma = 5e6
    Mw = catParent.data['Mag']
    M0 = 10 ** (1.5 * (Mw + 6.03))
    l = 0.001*(7 / 16 * M0 / dsigma) ** (1 / 3)
    ## hd
    HD = dNND['aHD']
    ## Rl
    aRl = HD / l
    print("haversine distance:",HD)
    print("rupture dimension",l)
    ###histogram
    aBins = np.arange(0, 30, dPar['eta_binsize'])
    aHist, aBins = np.histogram(aRl, aBins)
    aHist, aBins = np.array(list(zip(*zip(aHist, aBins))))  # cut to same length
    aRate_bin, aRate = seis_utils.eqRate(aRl,k)
    aHist /= dPar['eta_binsize']
    # =================================4==============================================
    #                          plot histogram
    # ================================================================================
    #myplot.plot_rate_pd(aBins, aHist, aRate_bin, aRate,curr_Mc,'a_hs',dPar['eta_binsize'])
    myplot.plot_loglog(aBins,aHist,aRate_bin,aRate,curr_Mc,'b_hs',k)




