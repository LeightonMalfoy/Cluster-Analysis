'''
plot the R/L or Haversine Distance distribution of aftershock sequences of certain mainshoks
-- background events excluded
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
#=================================0==============================================
#                            dir, file, params
#================================================================================
dir_in = '%s/data'%( os.path.expanduser('.'))
file_in= 'hs_1981_2018_all.mat'

#file_b  = '%s_b_Mc_D.txt'%(fileIn.split('.')[0])
dPar  = {   'aMc'         :  np.array( [2.5,2.6,3.0,3.5,4.0]),#np.array([3.0, 4.0]),np.array( [2.0, 2.5, 3.0, 3.5]
            # fractal dimension and b for eq. (1)
            'D'           : 1.6, #1.6  TODO: - these values should be contrained independently
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
curr_Mc = dPar['aMc'][0]
klist = [30,30,30]#150
minklst = [50]#80, ,
namelist=['Baja']#,'Northridge','HectorMine','Landers'

for name,k in zip(namelist,klist):
    # load nearest neighbor distances
    NND_file = 'data/NND_%s.mat' % name
    dNND   = data_utils.loadmat(NND_file) #,  struct_as_record=True)

    #load Landers/Joshua Tree earthquake catalog
    cat_file='data/cat_%s.mat' % name
    eqCat.loadMatBin(cat_file)

    # ==================================2=============================================
    #                       compute space-time-magnitude distance, histogram
    # ================================================================================
    # select only the clustering event pairs
    sel_cl = np.log10(dNND['aNND']) <= -5 #-4.7
    # HD
    HD = dNND['aHD'][sel_cl]
    print('catalog size: ', len(HD))
    #print('*************************', HD)
    ## l
    #l = 10 ** (-3.22 + 0.69 * eqCat.data['Mag'][0]) / 2 # unit: given by wells and coppersman
    #l=10**(0.44*eqCat.data['Mag'][0]-2) # given by AGU paper
    L = 1.99*10 ** (0.49 * eqCat.data['Mag'][0]-3)# calc by myself according to AGU
    ## Rl
    aRl = HD / L
    #======================================
    # histogram
    aBins = np.arange(0, 15, dPar['eta_binsize']) # TODO: change the bins limit
    aHist, aBins = np.histogram(aRl, aBins)
    aHist, aBins = np.array(list(zip(*zip(aHist,aBins))))  # cut to same length
    # correct for binsize
    aHist /= dPar['eta_binsize']
    aHist /= len(aRl)
    # ======================================
    # linear rate
    aRate_bin, aRate = seis_utils.eqlgRate(np.log10(aRl[aRl>0]),10)
    aRate /= len(aRl[aRl>0])
    # ======================================
    # cumulative distribution function
    aCum = np.cumsum(aHist)
    # ==================================3=============================================
    #                       get the mainshock infomation
    # ================================================================================
    Mag = eqCat.data['Mag'][0]
    YR = int(eqCat.data['Time'][0])
    # =================================4==============================================
    #                          NND dist
    # ================================================================================
    ###histogram
    aNND = np.log10(dNND['aNND'])
    aBins_NND = np.arange(-13, 1, dPar['nnd_binsize'])
    aHist_NND, aBins_NND = np.histogram(aNND, aBins_NND)
    aHist_NND, aBins_NND = np.array(list(zip(*zip(aHist_NND, aBins_NND))))  # cut to same length
    # correct for binsize
    aHist_NND /= dPar['nnd_binsize']
    # aHist /= len(aRl)
    aHist_NND /= len(aNND)
    # =================================5==============================================
    #                         haversine dist
    # ================================================================================
    aHD = dNND['aHD'][sel_cl]
    aBins_HD = np.arange(0, 500, dPar['hd_binsize'])
    aHist_HD, aBins_HD = np.histogram(aHD, aBins_HD)
    aHist_HD, aBins_HD = np.array(list(zip(*zip(aHist_HD, aBins_HD))))  # cut to same length
    # correct for binsize
    aHist_HD /= dPar['hd_binsize']
    # aHist /= len(aRl)
    aHist_HD /= len(aHD)
    aRateBin_HD, aRate_HD = seis_utils.eqlgRate(np.log10(aHD),k)
    aRate_HD /= len(aHD)
    '''
    # try: calc rate for lograthmic equally separated dot 
    aRateBin_HD = np.arange(np.amin(np.log10(aHD)),np.amax(np.log10(aHD)),0.005)
    aRate_HD, aRateBin_HD = np.histogram(np.log10(aHD),aRateBin_HD)
    aRateBin_HD = 10**aRateBin_HD
    binsize = aRateBin_HD[1:]-aRateBin_HD[:len(aRateBin_HD)-1]
    aRate_HD, aRateBin_HD = np.array(list(zip(*zip(aRate_HD, aRateBin_HD))))
    aRate_HD /= binsize
    aRate_HD /= len(aHD)
    aRateBin_HD = aRateBin_HD[aRate_HD>0]
    aRate_HD = aRate_HD[aRate_HD>0]
    '''
    # =================================6==============================================
    #                         cumulative distr
    # ================================================================================
    '''
    aHD = dNND['aHD'][sel_cl]
    aBins_HD = np.arange(0, 500, dPar['hd_binsize'])
    aHist_HD, aBins_HD = np.histogram(aHD, aBins_HD)
    aHist_HD, aBins_HD = np.array(list(zip(*zip(aHist_HD, aBins_HD))))  # cut to same length
    # correct for binsize
    aHist_HD /= dPar['hd_binsize']
    # aHist /= len(aRl)
    aHist_HD /= len(aHD)
    aRateBin_HD, aRate_HD = seis_utils.eqRate(aHD,k,dividing=1*L)#,minK=70
    aRate_HD /= len(aHD)
    # ======================================
    # cumulative distribution function
    aCum_HD = np.cumsum(aHist_HD)
    '''
    # =================================7==============================================
    #                          plot histogram
    # ================================================================================
    #myplot.plot_rate_pd(aBins, aHist, curr_Mc, 'h_RL1_%s' % name, dPar['eta_binsize'],aRate_bin=aRate_bin,aRate=aRate)
    #myplot.plot_rate_pd(aBins_NND, aHist_NND, curr_Mc, 'h_NND_%s' % name, dPar['nnd_binsize'])
    #myplot.plot_loglog(aBins, aHist, aRate_bin, aRate, curr_Mc, 'h_lgRL1_%s' % name, dPar['eta_binsize'])
    #myplot.plot_rate_pd(aBins_HD, aHist_HD, curr_Mc, name, dPar['hd_binsize'],Mag = Mag,YR=YR,L=L,plotkind='HD',aRate_bin=aRateBin_HD,aRate=aRate_HD)
    myplot.plot_loglog(aBins_HD, aHist_HD, curr_Mc, name, Mag = Mag, YR=YR, L=L, plotkind='lgHD4', aRate_bin=aRateBin_HD, aRate=aRate_HD,q=0.35,d=1.2,gamma=0.6,n=-1.35,c=1)#
    #myplot.plot_cumlog(aBins_HD, aHist_HD, aCum_HD, curr_Mc, name, Mag = Mag, YR=YR,  plotkind='HDcum')#

'''
for curr_Mc,k in zip(dPar['aMc'],dPar['k']):
    # cut below current completeness
    eqCatMc.copy(eqCat)
    eqCatMc.selectEvents(curr_Mc, None, 'Mag')


    # load nearest neighbor distances
    NND_file = 'data/%s_NND_Mc_%.1f_HD.mat'%( file_in.split('.')[0], curr_Mc)
    dNND   = data_utils.loadmat(NND_file) #,  struct_as_record=True)
    # ==================================2=============================================
    #                       compute space-time-magnitude distance, histogram
    # ================================================================================
    # only plot for clusterred events
    dCluster = {}
    eqCatCluster = EqCat()
    eqCatCluster.data['Mag'] = np.array([])
    dNND['aNND'] = np.log10(dNND['aNND'])
    catChild.copy(eqCatMc)
    catParent.copy(eqCatMc)
    # catChild, catPar = create_parent_child_cat( projCat, dNND)
    catChild.selEventsFromID(dNND['aEqID_c'], repeats=True)
    catParent.selEventsFromID(dNND['aEqID_p'], repeats=True)
    print('size of parent catalog', catParent.size(), 'size of offspring cat', catChild.size())
    # HD
    HD = np.array([])
    for iEv in list(range(catParent.size())):
        if dNND['aNND'][iEv] < dPar['eta_0']:  # cluster events
            eqCatCluster.data['Mag']=np.append(eqCatCluster.data['Mag'],catParent.data['Mag'][iEv])
            HD=np.append(HD, dNND['aHD'][iEv])
    print('current catalog size: ', len(HD))
    #print('*************************', HD)
    ## l
    sigma = 5e6
    Mw = eqCatCluster.data['Mag']
    M0 = 10 ** (1.5 * (Mw + 6.03))
    l = 0.001*(7 / 16 * M0 / sigma) ** (1 / 3)  # meter -> kilomiter
    ## Rl
    aRl = HD / l
    ###histogram
    aBins = np.arange(0, 20, dPar['eta_binsize'])
    aHist, aBins = np.histogram(aRl, aBins)
    aHist, aBins = np.array(list(zip(*zip(aHist, aBins))))  # cut to same length
    # correct for binsize
    aHist /= dPar['eta_binsize']
    aRate_bin, aRate = seis_utils.eqRate(aRl,k)
    # =================================4==============================================
    #                          plot histogram
    # ================================================================================
    #myplot.plot_rate_pd(aBins, aHist, aRate_bin, aRate, curr_Mc, 'c_cl', dPar['eta_binsize'])
    myplot.plot_loglog(aBins, aHist, aRate_bin, aRate, curr_Mc, 'd_cl', dPar['eta_binsize'])
'''