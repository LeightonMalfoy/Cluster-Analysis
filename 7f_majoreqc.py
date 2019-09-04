'''
Created on August 16, 2016

- compute inter-event times and distance and normalize by magnitude

@author: tgoebel
'''
import matplotlib as mpl

mpl.use('Agg')  # uncomment for interactive plotting

import sys
sys.path.append('/auto/home/lhuang/PycharmProjects/clustering-analysis-master')
import numpy as np
import matplotlib.pyplot as plt
from src.EqCat import *
from selector.getBigEQ import getMajor
from selector.getLanders import getLanders
from src.plot_major import *
import scipy.io

eqCat = EqCat()  # original catalog
eqCatMc = EqCat()  # this catalog wil be modfied with each Mc iteration
catChild = EqCat()
catParent = EqCat()
# =================================1==============================================
#                            dir, file, params
# ================================================================================
map_filename = 'plots/'
fault_file = 'data/CA_faults.all.dat'
# file_b  = '%s_b_Mc_D.txt'%(fileIn.split('.')[0])
dPar = {'aMc': np.array([2.5]),  # np.array([3.0, 4.0]),np.array( [2.0, 2.5, 3.0, 3.5]
        # fractal dimension and b for eq. (1)
        'D': 1.6,  # 1.6  TODO: - these values should be contrained independently
        'b': 1.0,
        # =================histogram==========
        'HD_binsize': 5,
        'eta_binsize': .05,
        }

for curr_Mc in dPar['aMc']:
    # =================================2==============================================
    #                            load data, select events
    # ================================================================================
    aCluster,dNND = getMajor(curr_Mc)

    for i in list(range(len(aCluster))): #nf len(aCluster)
        #i=1  # only for Landers
        cluster = aCluster[i]
        cluster.saveMatBin('../data/majorEQ_%d.mat' % i)
        # ====================================3===========================================
        #                       plot HD histogram
        # ================================================================================
        sel_c = cluster.data['sel_c'] # offspring
        aHD = dNND['aHD'][sel_c]
        ####### get NND file for the cluster
        NND_file = '../data/NND_majorEQ_%d.mat' % i
        nNND = {}
        nNND['aEqID_p'] = dNND['aEqID_p'][sel_c]
        nNND['aEqID_c'] = dNND['aEqID_c'][sel_c]
        nNND['aNND'] = dNND['aNND'][sel_c]
        nNND['aHD'] = aHD
        scipy.io.savemat(NND_file, nNND, do_compression=True)
        #######
        sel_cluster = np.log10(dNND['aNND'][sel_c]) <= -5
        aHD1 = dNND['aHD'][sel_c][sel_cluster]
        HD_file = '../plots/g1_HD_his_%d.png' % i
        majorEQ = EqCat()
        majorEQ.data['Time']=cluster.data['Time'][0]
        majorEQ.data['Mag']=cluster.data['Mag'][0]
        #plot_HD(HD_file, aHD, majorEQ, dPar['HD_binsize'], curr_Mc)
        # ====================================3'===========================================
        #plot_HD1(HD_file, aHD,aHD1, majorEQ, dPar['HD_binsize'], curr_Mc)
        # ====================================4===========================================
        #                       plot NND histogram
        # ================================================================================
        aEta = np.log10(dNND['aNND'][sel_c])
        eta_file = '../plots/g_NND_his_%d.png' % i

        #plot_eta(eta_file, aEta, majorEQ, dPar['eta_binsize'], curr_Mc)
        # ====================================5===========================================
        #                       plot lgR-lgτ distribution
        # ================================================================================
        """
        aEta1 = np.array([])
        aR = np.array([])
        aTau = np.array([])
        # calculate the lgR, lgτ
        mstime=cluster.data['Time'][0]
        msmag=cluster.data['Mag'][0]
        msx=cluster.data['X'][0]
        msy=cluster.data['Y'][0]
        print("main shock magnitude:%.1f"%msmag)
        for x, y, time in zip(cluster.data['X'][1:], cluster.data['Y'][1:], cluster.data['Time'][1:]):
            t = (time - mstime) * 10 ** (-0.5*dPar['b'] * msmag)
            r = ((x - msx) ** 2 + (y - msy) ** 2) ** (dPar['D'] / 2) * 10 ** (-0.5*dPar['b'] * msmag)
            eta = r * t
            aEta1 = np.append(aEta1, np.log10(eta))
            aTau = np.append(aTau, np.log10(t))
            aR = np.append(aR, np.log10(r))
        dPar1 = {'binx' : 0.1, 'biny' : .1,# used for density and gaussian smoothing
                        'sigma'   : None, #.2, #'silverman', # Gaussian smoothing bandwidth, default = n**(-1./(d+4)),
                    #=================plotting==============
                    'eta_0'       : -5.0,
                    #'xmin' : -13, 'xmax' : 0,
                    'Tmin' :  -8, 'Tmax' : 0,
                    'Rmin' :  -5, 'Rmax' : 3,
                    'cmap' : plt.cm.RdYlGn_r, }#aofRT brg
        rtau_file= '../plots/g_tau_his_%d.png' % i
        ____,at1,ar1,___=getLanders()

        '''dfgfdah'''
        a_Tbin = np.arange(dPar1['Tmin'], dPar1['Tmax'] + 2 * dPar1['binx'], dPar1['binx'])
        a_Rbin = np.arange(dPar1['Rmin'], dPar1['Rmax'] + 2 * dPar1['biny'], dPar1['biny'])
        XX, YY, ZZ = density_2D(aTau, aR, a_Tbin, a_Rbin, sigma=dPar1['sigma'])
        plt.figure(1, figsize=(8, 10))
        ax = plt.subplot(111)
        normZZ = ZZ * (dPar1['binx'] * dPar1['biny'] * len(aTau))
        plot1 = ax.pcolormesh(XX, YY, normZZ, cmap=dPar1['cmap'])
        cbar = plt.colorbar(plot1, orientation='horizontal', shrink=.5, aspect=20, )
        # plot eta_0 to divide clustered and background mode
        ax.plot([dPar1['Tmin'], dPar1['Tmax']], -np.array([dPar1['Tmin'], dPar1['Tmax']]) + dPar1['eta_0'], '-', lw=1.5,
                color='w')
        ax.plot([dPar1['Tmin'], dPar1['Tmax']], -np.array([dPar1['Tmin'], dPar1['Tmax']]) + dPar1['eta_0'], '--', lw=1.5,
                color='.5')
        # plot eta_0 to divide clustered and background mode
        x = np.array([dPar1['Tmin'], dPar1['Tmax']])
        y = 1.6*x+14.952
        ax.plot(x, y, '--', lw=1.5,
                color='.5')
        cbar.set_label('Number of Event Pairs', labelpad=-40)
        #ax.plot(at1,ar1,'bo',color=None)
        ax.set_xlabel('Rescaled Time')
        ax.set_ylabel('Rescaled Distance')
        ax.set_ylim([-5.0,3.0])
        ax.set_title("Year:%d Mag:%.1f/nNearest Neighbor Pairs in R-T" % (majorEQ.data['Time'], majorEQ.data['Mag']))
        plt.savefig('../plots/h_Landers_RTau.png', dpi=500)
        plt.clf()
        ###
        """


        #plot_r_tau_density(rtau_file, aTau, aR, dParaofRT,majorEQ)