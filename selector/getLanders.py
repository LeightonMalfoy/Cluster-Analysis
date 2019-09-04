import sys
sys.path.append('/auto/home/lhuang/PycharmProjects/clustering-analysis-master')
import matplotlib as mpl
# mpl.use( 'Agg')
import matplotlib.pyplot as plt

import numpy as np
import os

# ------------------------------my modules--------------------------------------
import src.data_utils as data_utils
import src.clustering as clustering
from src.EqCat import *
import src.seis_utils as seis_utils



def getLanders():
    eqCat = EqCat()  # original catalog
    eqCatMc = EqCat()  # this catalog will be modified with each Mc iteration
    catLanders= EqCat()
    catCoso = EqCat()

    # =================================1==============================================
    #                            dir, file, params
    # ================================================================================
    data_dir = '../data'
    plot_dir = '../plots'
    file_in = 'hs_1981_2011_all.mat'

    dPar = {'a_Mc': np.array([3.0, 4.0]),  # np.array( [2.0, 2.5, 3.0, 3.5]),
            # separate clustered and background
            'eta_0': -5.0,
            'testPlot': True,

            'D': 1.6,  # 1.6  TODO: - these values should be contrained independently
            'b': 1.0,
            }

    curr_Mc = dPar['a_Mc'][0]
    # =================================2==============================================
    #                            load data, select events
    # ================================================================================
    eqCat.loadMatBin(os.path.join(data_dir, file_in))
    print('total no. of events', eqCat.size())
    eqCat.selectEvents(curr_Mc, None, 'Mag')
    #eqCat.toCart_coordinates()
    # eqCat.selector( tmin, tmax, 'Time')
    print('no. of events after initial selection', eqCat.size())

    catLanders.copy(eqCat)
    catCoso.copy(eqCat)

    catLanders.selectEvents(7.0,None,'Mag')
    catLanders.selectEvents(34,35,'Lat')
    catLanders.selectEvents(-117,-116,'Lon')
    catLanders.selectEvents(1992,1993,'Time')
    #print("===========Landers Info============\n",catLanders.data)

    catCoso.selectEvents(catLanders.data['Time'][0],catLanders.data['Time'][0]+0.1,'Time')
    catCoso.selectEvents(-118.5,-117,'Lon')
    catCoso.selectEvents(35.5,36.5,'Lat')
    #print("===========Coso Info===============\nLon\tLat\tMag\tTime\t")
    #for lon,lat,mag,time in zip(catCoso.data['Lon'],catCoso.data['Lat'],catCoso.data['Mag'],catCoso.data['Time']):
    #    print("%.3f\t%.3f\t%.2f\t%.8f\t"%(lon,lat,mag,time))

    aEta = np.array([])
    aT = np.array([])
    aR = np.array([])

    catAll=EqCat()
    catAll.merge(catLanders,catCoso)
    catAll.toCart_coordinates()

    #print(catAll.data)

    for x,y,time in zip(catAll.data['X'][1:],catAll.data['Y'][1:],catAll.data['Time'][1:]):
        x=catAll.data['X'][0]
        y=catAll.data['Y'][0]+3
        t = (time-catAll.data['Time'][0])*10**(-dPar['b']*catAll.data['Mag'][0]/2)
        print("distance:",((x-catAll.data['X'][0])**2+(y-catAll.data['Y'][0])**2)**0.5)
        r = ((x-catAll.data['X'][0])**2+(y-catAll.data['Y'][0])**2)**(dPar['D']/2)*10**(-dPar['b']*catAll.data['Mag'][0]/2)
        eta=r*t
        aEta=np.append(aEta, np.log10(eta))
        aT=np.append(aT, np.log10(t))
        aR=np.append(aR, np.log10(r))
    print("===========Neareast Neighbor Distance===============\nr\tt\teta\t")
    for r, t, eta in zip(aR,aT , aEta):#(catCoso.data['Time']-catLanders.data['Time'][0])*365.25*24*3600
        print("%f\t%f\t%f" % (r, t, eta))

    return catAll,aT,aR,aEta

getLanders()
