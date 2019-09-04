
'''
- map the distribution of aftershocks of certain major earthquakes
    colors show
    1. the Δ lg(NND)(nearest neighbor distance) between aftershock and its parent
    or,
    2. the interval time between aftershock and its parent
- input file:
    1. catalog of certain earthquakes,

@author: Litong Huang, Peking University
'''

import matplotlib as mpl
mpl.use('Agg')  # uncomment for interactive plotting

import os
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
#import matplotlib.colors as co
import src.data_utils as data_utils
#import src.seis_utils as seis_utils
from src.EqCat import *
#from selector.getBigEQ import getMajor

eqCat = EqCat()  # original catalog
eqCatMc = EqCat()  # this catalog wil be modfied with each Mc iteration
catChild = EqCat()
catParent = EqCat()
# =================================1==============================================
#                            dir, file, params
# ================================================================================
fault_file = 'data/CA_faults.all.dat'
# file_b  = '%s_b_Mc_D.txt'%(fileIn.split('.')[0])
plottype = 'time' # TODO: change here for different plot purpose. or 'eta'
dPar = {'aMc': np.array([2.5]),  # np.array([3.0, 4.0]),np.array( [2.0, 2.5, 3.0, 3.5]
        # fractal dimension and b for eq. (1)
        'D': 1.6,  # 1.6  TODO: - these values should be contrained independently
        'b': 1.0,
        # =================map plotting==========
        'minlon': -124, 'minlat': 30,
        'maxlon': -114, 'maxlat': 39,
        'minmag': 2.5,
        'R_Lmin': 25,  # to select event-pairs in the tail
        }

for curr_Mc in dPar['aMc']:
    # =================================2==============================================
    #                            load data, select events
    # ================================================================================
    #aCluster,dNND = getMajor(curr_Mc)
    """
    for i in list(range(len(aCluster))):#
        i=4
        cluster = aCluster[i]
        # ================================================================================
        #                       calculate delta NND for each pair
        # ================================================================================
        sel_c = cluster.data['sel_c']
        adNND = np.log10(dNND['aNND'][sel_c]) + 5
        #print('min Δ lg(NND):',min(adNND),'max Δ lg(NND):' ,max(adNND))
        ### select only the cluster events # TODO: example - dNND['aNND'][sel_c][sel_cluster], adNND[sel_cluster], adTime[sel_cluister]
        sel_cluster = adNND <= 0 """

    namelist = ['Baja']#,'HectorMine','Landers'
    for name in namelist:
        # load nearest neighbor distances
        NND_file = 'data/NND_%s.mat' % name
        # NND_file = 'data/%s_NND_Mc_%.1f_HD.mat'%( file_in.split('.')[0], curr_Mc)
        dNND = data_utils.loadmat(NND_file)  # ,  struct_as_record=True)

        # load Landers/Joshua Tree earthquake catalog
        cat_file = 'data/cat_%s.mat' % name
        eqCat.loadMatBin(cat_file)
        # ================================================================================
        #                       calculate interval time for each pair
        # ================================================================================
        adTime = eqCat.data['Time'][1:]-eqCat.data['Time'][0]
        aNND = dNND['aNND']
        sel_cluster = np.log10(dNND['aNND']) <= -5
        aClTime = adTime[sel_cluster]
        adNND = aNND[sel_cluster] + 5
        cluster = EqCat()
        cluster.copy(eqCat)
        print('all clustering aftershocks:',len(aClTime))
        if len(aClTime)<=200:
            print("***************************************************")
            print("%s, insufficient offsprings." % name)
        else:
            print("***************************************************")
            print("****%s, sufficient offsprings." % name)
            sel_time = aClTime <= 3/365.25
            sel_time1 = aClTime > 3/365.25
            aClTime = aClTime[sel_time]
            aClTime *= 365.25*24 # --> hr
            print('min time interval:',min(aClTime),'max time interval:' ,max(aClTime))
            print('earthquakes within 72h:',len(aClTime))
            # ================================================================================
            #                 plot map
            # ================================================================================
            fig = plt.figure()
            ax = plt.subplot()

            # set basemap
            mymap = Basemap(llcrnrlon=dPar['minlon'], llcrnrlat=dPar['minlat'],
                            urcrnrlon=dPar['maxlon'], urcrnrlat=dPar['maxlat'],
                            resolution='f')
            mymap.drawcoastlines(linewidth=0.25)
            mymap.drawcountries(linewidth=0.25)
            mymap.drawstates(linewidth=0.25)
            mymap.fillcontinents(color='none', lake_color='aqua')

            # plot faults
            data = np.loadtxt(fault_file)
            x = []
            y = []
            for plot in data:
                if plot[0] < 0:
                    x.append(plot[0])
                    y.append(plot[1])
                else:
                    mymap.plot(x, y, marker='None', color='black', linewidth=.2)
                    x = []
                    y = []

            # labels = [left,right,top,bottom]
            parallels = np.arange(dPar['minlat'], dPar['maxlat'], 2)
            mymap.drawparallels(parallels, labels=[True, False, True, False])
            meridians = np.arange(dPar['minlon'], dPar['maxlon'], 2)
            mymap.drawmeridians(meridians, labels=[True, False, False, True])

            cmap = plt.cm.RdYlGn

            X_Child, Y_Child = mymap(cluster.data['Lon'][1:][sel_cluster][sel_time],
                                     cluster.data['Lat'][1:][sel_cluster][sel_time])
            X_oChild, Y_oChild = mymap(cluster.data['Lon'][1:][sel_cluster][sel_time1],
                                       cluster.data['Lat'][1:][sel_cluster][sel_time1])
            X_Parent, Y_Parent = mymap(cluster.data['Lon'][0], cluster.data['Lat'][0])

            # plot events
            if plottype == 'time':
                ### color represent time interval
                eqChildrenOther = mymap.scatter(X_oChild, Y_oChild,
                                                .4 * cluster.data['Mag'][1:][sel_cluster][sel_time1],
                                                marker='o', \
                                                c='.5', alpha=1,
                                                label='Eq after 72 hrs: %d' % len(X_oChild))  # c=times,
                eqChildren = mymap.scatter(X_Child, Y_Child, .4 * cluster.data['Mag'][1:][sel_cluster][sel_time],
                                           marker='o', \
                                           c=aClTime, vmin=min(aClTime), vmax=max(aClTime), cmap=cmap, alpha=1,
                                           label='Eq within 72 hrs: %d' % len(X_Child))  # c=times,
            if plottype == 'eta':
                ### color represent Δ lg(NND)
                eqChildrenOther = mymap.scatter(X_oChild, Y_oChild,
                                                .4 * cluster.data['Mag'][1:][sel_cluster][sel_time1],
                                                marker='o', \
                                                c='.5', alpha=1,
                                                label='Eq after 72 hrs: %d' % len(X_oChild))  # c=times,
                eqChildren = mymap.scatter(X_Child, Y_Child, .4 * cluster.data['Mag'][1:][sel_cluster][sel_time],
                                           marker='o', \
                                           c=adNND, vmin=min(adNND), vmax=-min(adNND), cmap=cmap, alpha=1,
                                           label='Eq within 72 hrs: %d' % len(X_Child))  # c=times,

            plt.legend()
            # eqChildren = mymap.scatter(np.flip(X_Child,axis=0), np.flip(Y_Child,axis=0), np.flip(.5*cluster.data['Mag'][1:][sel_cluster], axis=0), marker='o', c=np.flip(adNND[sel_cluster],axis=0),vmin=min(adNND),vmax=0,cmap=cmap, alpha=1)  # c=times,
            eqParents1 = mymap.scatter(X_Parent, Y_Parent, 2*cluster.data['Mag'][0], marker='*', c='black', alpha=1)  # c=times,
            eqParents2 = mymap.scatter(X_Parent, Y_Parent, 1*cluster.data['Mag'][0], marker='*', c='red', alpha=1)  # c=times,

            if plottype == 'time':
                #### colorbar >> colors represent interval time
                cbar=plt.colorbar(eqChildren)
                cbar.set_label("interval time/hr")
                ticklist = np.arange(int(min(aClTime))+10,int(max(aClTime)), 10)
                #ticklist = np.arange(int(min(adNND)),-int(min(adNND)),1)
                cbar.set_ticks(ticklist)#
                cbar.set_ticklabels(ticklist)

            if plottype == 'eta':
                #### colorbar >> colors represent Δ lg(NND)
                cbar=plt.colorbar(eqChildren)
                cbar.set_label(" Δ lg(NND)")
                ticklist = np.arange(int(min(adNND)),0, 1)
                #ticklist = np.arange(int(min(adNND)),-int(min(adNND)),1)
                cbar.set_ticks(ticklist)
                cbar.set_ticklabels(ticklist)

            ax.set_title("%s,Year:%d, Mag:%.1f" % (name, np.floor(cluster.data['Time'][0]),cluster.data['Mag'][0]))

            plt.show()
            fig.savefig('plots/h_map_%s_%s.png'% (name,plottype), dpi=500) # it -- interval time
            plt.clf()
