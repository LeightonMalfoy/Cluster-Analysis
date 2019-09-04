'''
Created on August 16, 2016

- compute inter-event times and distance and normalize by magnitude

@author: tgoebel
'''
import matplotlib as mpl

mpl.use('Agg')  # uncomment for interactive plotting

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
# ------------------------------my modules--------------------------------------
import src.data_utils as data_utils
import src.clustering as clustering
from src.EqCat import *

eqCat = EqCat()  # original catalog
eqCatMc = EqCat()  # this catalog wil be modfied with each Mc iteration
catChild = EqCat()
catParent = EqCat()
# =================================1==============================================
#                            dir, file, params
# ================================================================================
dir_in = '%s/data' % (os.path.expanduser('.'))
file_in = 'hs_1981_2018_all.mat'

# file_b  = '%s_b_Mc_D.txt'%(fileIn.split('.')[0])
dPar = {'a_Mc': np.array([2.5, 3.0, 3.5, 4.0]),  # np.array([3.0, 4.0]),
        # fractal dimension and b for eq. (1)
        'D': 1.6,  # 1.6  TODO: - these values should be constrained independently
        'b': 1.0,
        'M_pt': 3.5, # threshold :only select event pairs with parent event larger than it

        # ===================smoothing parameters=============
        'binx': .1, 'biny': .1,  # used for density and gaussian smoothing
        'sigma': None,  # .2, #'silverman', # Gaussian smoothing bandwidth, default = n**(-1./(d+4)),
        # =================plotting==============
        'eta_0': -5.0,
        # 'xmin' : -13, 'xmax' : 0,
        'Tmin': -8, 'Tmax': 0,
        'Rmin': -5, 'Rmax': 3,
        'cmap': plt.cm.RdYlGn_r,
        'R_L': 25,
        'R_Lmax':5,}

# =================================2==============================================
#                            load data, select events
# ================================================================================
eqCat.loadMatBin(os.path.join(dir_in, file_in))
print('total no. of events', eqCat.size())
eqCat.selectEvents(dPar['a_Mc'][0], None, 'Mag')
# eqCat.selector( tmin, tmax, 'Time')
print('no. of events after initial selection', eqCat.size())

# two ways to do the distance comp: 1 project into equal distance azimuthal , comp Cartersian distance in 3D
#                                   2 get surface distance from lon, lat (haversine), use pythagoras to include depth
eqCat.toCart_coordinates(projection='aeqd')

for i in range(dPar['a_Mc'].shape[0]):
    f_Mc = dPar['a_Mc'][i]

    # cut below current completeness
    eqCatMc.copy(eqCat)
    eqCatMc.selectEvents(f_Mc, None, 'Mag')
    print('current catalog size: ', eqCatMc.size())

    # load nearest neighbor distances
    NND_file = 'data/%s_NND_Mc_%.1f_HD.mat' % (file_in.split('.')[0], f_Mc)
    dNND = data_utils.loadmat(NND_file)  # ,  struct_as_record=True)
    dNND['aNND'] = np.log10(dNND['aNND'])
    # ==================================3.0 ==========================================
    #                    select event by Mc, parent MC
    # ================================================================================
    catChild.copy(eqCatMc)
    catParent.copy(eqCatMc)
    # catChild, catPar = create_parent_child_cat( projCat, dNND)
    catChild.selEventsFromID(dNND['aEqID_c'], repeats=True)
    catParent.selEventsFromID(dNND['aEqID_p'], repeats=True)
    print('before::: size of parent catalog', catParent.size(), 'size of offspring cat', catChild.size())
    # select event pairs with parent event larger than M_pt
    sel = catParent.data['Mag'] >= dPar['M_pt']
    catChild.selEventsFromID(dNND['aEqID_c'][sel], repeats=True)
    catParent.selEventsFromID(dNND['aEqID_p'][sel], repeats=True)
    print('after::: size of parent catalog', catParent.size(), 'size of offspring cat', catChild.size())

    # ==================================3=============================================
    #                       compute re-scaled interevent times and distances
    # ================================================================================

    # note that dictionary dPar here has to include 'b','D' and 'Mc'
    a_R, a_T = clustering.rescaled_t_r(catChild, catParent, {'b': dPar['b'], 'D': dPar['D'], 'Mc': f_Mc},
                                       correct_co_located=True)
    # ==================================4==============================================================
    #                       T-R density plots
    # =================================================================================================
    a_Tbin = np.arange(dPar['Tmin'], dPar['Tmax'] + 2 * dPar['binx'], dPar['binx'])
    a_Rbin = np.arange(dPar['Rmin'], dPar['Rmax'] + 2 * dPar['biny'], dPar['biny'])
    XX, YY, ZZ = data_utils.density_2D(np.log10(a_T), np.log10(a_R), a_Tbin, a_Rbin, sigma=dPar['sigma'])

    plt.figure(1, figsize=(8, 16))
    ax = plt.subplot(211)
    ax.set_title('Nearest Neighbor Pairs in R-T')
    # ------------------------------------------------------------------------------
    normZZ = ZZ * (dPar['binx'] * dPar['biny'] * eqCatMc.size())
    plot1 = ax.pcolormesh(XX, YY, normZZ, cmap=dPar['cmap'])
    cbar = plt.colorbar(plot1, orientation='horizontal', shrink=.5, aspect=20, )
    # ax.plot(  np.log10( a_T), np.log10( a_R), 'wo', ms = 1.5, alpha = .2)
    # plot eta_0 to divide clustered and background mode
    ax.plot([dPar['Tmin'], dPar['Tmax']], -np.array([dPar['Tmin'], dPar['Tmax']]) + dPar['eta_0'], '-', lw=1.5,
            color='w')
    ax.plot([dPar['Tmin'], dPar['Tmax']], -np.array([dPar['Tmin'], dPar['Tmax']]) + dPar['eta_0'], '--', lw=1.5,
            color='.5')

    # -----------------------labels and legends-------------------------------------------------------
    # cbar.set_label( 'Event Pair Density [#ev./dRdT]')
    cbar.set_label('Number of Event Pairs', labelpad=-40)
    ax.set_xlabel('Rescaled Time')
    ax.set_ylabel('Rescaled Distance')
    ax.set_xlim(dPar['Tmin'], dPar['Tmax'])
    ax.set_ylim(dPar['Rmin'], dPar['Rmax'])

    # ==================================5==============================================================
    #                       Large-R/L plots
    # =================================================================================================
    # compute R/L
    # HD
    HD = np.array([])
    catCluster = EqCat()
    catCluster.data['Mag'] = np.array([])
    for iEv in list(range(catChild.size())):
        #if dNND['aNND'][iEv] < dPar['eta_0']:  # cluster events
        #    catCluster.data['Mag'] = np.append(catCluster.data['Mag'], eqCatMc.data['Mag'][iEv + 1])
        #    HD = np.append(HD, dNND['aHD'][iEv])
        #catChild.data['Mag'] = np.append(catChild.data['Mag'], eqCatMc.data['Mag'][iEv + 1])
        HD = np.append(HD, dNND['aHD'][iEv])
    print('current catalog size: ', len(HD))
    # print('*************************', HD)
    ## l
    sigma = 5e6
    Mw = catParent.data['Mag']
    M0 = 10 ** (1.5 * (Mw + 6.03))
    l = 0.001 * (7 / 16 * M0 / sigma) ** (1 / 3)  # meter -> kilomiter
    ## Rl
    aRl = HD / l
    ###
    sel_TR = aRl >= dPar['R_L']
    # compute the density of Large-R/L events
    a_T1 = a_T[sel_TR]
    a_R1 = a_R[sel_TR]
    XX1, YY1, ZZ1 = data_utils.density_2D(np.log10(a_T1), np.log10(a_R1), a_Tbin, a_Rbin, sigma=dPar['sigma'])
    ay = plt.subplot(212)
    ay.set_title('R/L>=%.1f event pairs in R-T' % dPar['R_L'])
    # ------------------------------------------------------------------------------
    normZZ1 = ZZ1 * (dPar['binx'] * dPar['biny'] * eqCatMc.size())
    plot2 = ay.pcolormesh(XX1, YY1, normZZ1, cmap=dPar['cmap'])
    cbar2 = plt.colorbar(plot2, orientation='horizontal', shrink=.5, aspect=20, )
    cbar.set_label('Number of Event Pairs with parent bigger than %.1f' % dPar['M_pt'], labelpad=-40)
    ay.set_xlabel('Rescaled Time')
    ay.set_ylabel('Rescaled Distance')
    ay.set_xlim(dPar['Tmin'], dPar['Tmax'])
    ay.set_ylim(dPar['Rmin'], dPar['Rmax'])
    ay.plot([dPar['Tmin'], dPar['Tmax']], -np.array([dPar['Tmin'], dPar['Tmax']]) + dPar['eta_0'], '-', lw=1.5,
            color='w')
    ay.plot([dPar['Tmin'], dPar['Tmax']], -np.array([dPar['Tmin'], dPar['Tmax']]) + dPar['eta_0'], '--', lw=1.5,
            color='.5')

    print('large HD event number:',len(a_R1))
    #ax.scatter(np.log10(a_T1), np.log10(a_R1), marker='s',color='black')
    # =================================5==============================================
    #                           save results
    # ================================================================================
    print('plot saved in: ', 'plots/T_R_%s_Mc_%.1f.png' % (file_in.split('.')[0], f_Mc))
    plt.figure(1)
    plt.savefig('plots/f_RL_%.1f_Mc_%.1f.png' % (dPar['R_L'], f_Mc))
    plt.show()
    plt.clf()