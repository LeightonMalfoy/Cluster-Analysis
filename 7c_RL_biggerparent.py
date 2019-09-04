'''
Created on August 16, 2016

- compute inter-event times and distance and normalize by magnitude

@author: tgoebel
'''
import matplotlib as mpl
mpl.use( 'Agg') # uncomment for interactive plotting

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
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
dPar  = {   'aMc'         :  np.array( [2.5]),#np.array([3.0, 4.0]),np.array( [2.0, 2.5, 3.0, 3.5]2.5,3.0,3.5,,  3.0,3.5
            # fractal dimension and b for eq. (1)
            'D'           : 1.6, #1.6  TODO: - these values should be contrained independently
            'b'           : 1.0,
            'M_pt': 6.5,  # threshold :only select event pairs with parent event larger than it
            #=================plotting==============
            'eta_binsize' :  .1                                                                                                                                                      ,
            'xmin' : 0, 'xmax' : 1.0,
            'eta_0': -5.0,
            #'k' :[int(300*(i+1)**(np.log10(1/6))) for i in list(range(10))],
            'k':[300], # to calculate the r/l rate  ,200,100
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

for curr_Mc,k in zip(dPar['aMc'],dPar['k']):
    # cut below current completeness
    eqCatMc.copy(eqCat)
    eqCatMc.selectEvents(curr_Mc, None, 'Mag')
    print( 'current catalog size: ',eqCatMc.size())

    # load nearest neighbor distances
    NND_file = 'data/%s_NND_Mc_%.1f_HD.mat'%( file_in.split('.')[0], curr_Mc)
    dNND   = data_utils.loadmat(NND_file) #,  struct_as_record=True)
    # ================================================================================
    #                    all event pairs
    # ================================================================================
    catChild.copy(eqCatMc)
    catParent.copy(eqCatMc)
    # catChild, catPar = create_parent_child_cat( projCat, dNND)
    catChild.selEventsFromID(dNND['aEqID_c'], repeats=True)
    catParent.selEventsFromID(dNND['aEqID_p'], repeats=True)
    print('before::: size of parent catalog', catParent.size(), 'size of offspring cat', catChild.size())

    dsigma = 5e6
    Mw = catParent.data['Mag']
    M0 = 10 ** (1.5 * (Mw + 6.03))
    l = 0.001 * (7 / 16 * M0 / dsigma) ** (1 / 3)
    ## hd
    HD = dNND['aHD']
    ## Rl
    aRl = HD / l
    #print("haversine distance:", HD)
    #print("rupture dimension", l)
    ###histogram
    aBins = np.arange(0, 30, dPar['eta_binsize'])
    aHist, aBins = np.histogram(aRl, aBins)
    aHist, aBins = np.array(list(zip(*zip(aHist, aBins))))  # cut to same length
    aRate_bin, aRate = seis_utils.eqRate(aRl, k)
    aHist /= dPar['eta_binsize']
    # ================================================================================
    #                       bigger parent event pairs
    # ================================================================================
    # select event pairs with parent event larger than M_pt
    sel = catParent.data['Mag'] >= dPar['M_pt']
    catChild.selEventsFromID(dNND['aEqID_c'][sel], repeats=True)
    catParent.selEventsFromID(dNND['aEqID_p'][sel], repeats=True)
    print('after::: size of parent catalog', catParent.size(), 'size of offspring cat', catChild.size())
    ## l
    Mw1 = catParent.data['Mag']
    M01 = 10 ** (1.5 * (Mw1 + 6.03))
    l1 = 0.001*(7 / 16 * M01 / dsigma) ** (1 / 3)
    ## hd
    HD1 = dNND['aHD'][sel]
    ## Rl
    aRl1 = HD1 / l1
    print("haversine distance:",HD)
    print("rupture dimension",l)
    ###histogram
    aBins1 = np.arange(0, 30, dPar['eta_binsize'])
    aHist1, aBins1 = np.histogram(aRl1, aBins1)
    aHist1, aBins1 = np.array(list(zip(*zip(aHist1, aBins1))))  # cut to same length
    aRate_bin1, aRate1 = seis_utils.eqRate(aRl1,k)
    aHist1 /= dPar['eta_binsize']
    aHist1 /= len(aRl1)
    aRate1 /= len(aRl1)
    # ================================================================================
    #                       plot - difference
    # ================================================================================
    '''#
    ax = plt.subplot()
    # ax.plot( vBin, vHist, 'ko')
    sel_Rate = [np.where(abs(aRate_bin-x) < 0.05)[0][0] for x in aRate_bin1]
    print(aRate_bin[sel_Rate]-aRate_bin1)
    ax.plot(aRate_bin1, aRate[sel_Rate]-aRate1, 'r-', label='Eq rate difference')
    ax.legend(loc='upper left')
    ax.set_xlabel('R/l')
    ax.set_ylabel('difference')
    ax.grid('on')
    plt.show()
    plotFile = 'plots/%s_Mc_%.1f.png' % ('a3', curr_Mc)
    print('save plot', plotFile)
    plt.savefig(plotFile)
    plt.clf()
    '''#
    # ================================================================================
    #                       plot -lglg
    # ================================================================================
    #'''
    ax = plt.subplot()
    # ax.plot( vBin, vHist, 'ko')
    ax.loglog(aBins1, aHist1, 'ko', mfc = 'none', mew = .2, label='Historical Rate')
    #ax.loglog(aRate_bin, aRate, 'r-', label='Earthquake rate - all pairs')
    ax.loglog(aRate_bin1, aRate1, 'b-', label='Earthquake rate - bigger parent than 6.5')
    ax.legend(loc='upper left')
    ax.set_xlabel('R/l')
    ax.set_ylabel('Probability Density')
    ax.grid('on')
    ax.set_title('Mc = %.1f,k=%d' % (curr_Mc,k))
    plt.show()
    plotFile = 'plots/%s.png' % ('h_major_lgRl')
    print('save plot', plotFile)
    plt.savefig(plotFile)
    plt.clf()
    #'''
    # ================================================================================
    #                       plot - distr
    # ================================================================================
    '''
    ax = plt.subplot()
    # ax.plot( vBin, vHist, 'ko')
    #ax.bar(aBins, aHist, width=0.9 * dPar['eta_binsize'], align='edge', color='.5')
    ax.bar(aBins1, aHist1, width=0.9 * dPar['eta_binsize'], align='edge', color='.5')
    #ax.plot(aRate_bin, aRate, 'b--', label='All event pairs, mc=2.5')
    ax.plot(aRate_bin1, aRate1, 'r-', label='Linear rate')
    ax.legend(loc='upper left')
    ax.set_xlabel('R/l')
    ax.set_ylabel('Probability Density')
    ax.set_xlim([0,30])
    ax.grid('on')
    plt.show()
    plotFile = 'plots/%s.png' % ('h_major_RL')
    print('save plot', plotFile)
    plt.savefig(plotFile)
    plt.clf()
    '''

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