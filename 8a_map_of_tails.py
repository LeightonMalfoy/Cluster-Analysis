'''
- map the distribution of event pairs in the R/L distribution tail

@author: Litong Huang
'''
import matplotlib as mpl
mpl.use( 'Agg') # uncomment for interactive plotting

import os
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.colors as co
import src.data_utils as data_utils
import src.seis_utils as seis_utils
from src.EqCat import *
from selector.getLanders import getLanders

eqCat   = EqCat() # original catalog
eqCatMc = EqCat() # this catalog wil be modfied with each Mc iteration
catChild=  EqCat()
catParent= EqCat()
#=================================1==============================================
#                            dir, file, params
#================================================================================
dir_in = '%s/data'%( os.path.expanduser('.'))
file_in= 'hs_1981_2018_all.mat'
map_filename = 'plots/'
fault_file = 'data/CA_faults.all.dat'
#file_b  = '%s_b_Mc_D.txt'%(fileIn.split('.')[0])
dPar  = {   'aMc'         :  np.array( [2.5,3.0,3.5,4.0]),#np.array([3.0, 4.0]),np.array( [2.0, 2.5, 3.0, 3.5]
            # fractal dimension and b for eq. (1)
            'D'           : 1.6, #1.6  TODO: - these values should be contrained independently
            'b'           : 1.0,
            'M_pt': 7.0,  # threshold :only select event pairs with parent event larger than it
            #=================plotting==============
            'eta_binsize' :  .1                                                                                                                                                      ,
            'xmin' : 0, 'xmax' : 1.0,
            'eta_0': -5.0,
            #'k' :[int(300*(i+1)**(np.log10(1/6))) for i in list(range(10))],
            'k':[600,400,200,50], # to calculate the r/l rate
            #=================map plotting==========
            'minlon':-123,'minlat':30,
            'maxlon':-113, 'maxlat':40,
            'minmag':2.5,
            'R_Lmin': 25, # to select event-pairs in the tail
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

#################################################
curr_Mc=dPar['aMc'][0]
k=dPar['k'][2]
# cut below current completeness
eqCatMc.copy(eqCat)
eqCatMc.selectEvents(curr_Mc, None, 'Mag')
print('current catalog size: ', eqCatMc.size())

# load nearest neighbor distances
NND_file = 'data/%s_NND_Mc_%.1f_HD.mat' % (file_in.split('.')[0], curr_Mc)
dNND = data_utils.loadmat(NND_file)  # ,  struct_as_record=True)
# ================================================================================
#                    all event pairs
# ================================================================================
catChild.copy(eqCatMc)
catParent.copy(eqCatMc)
# catChild, catPar = create_parent_child_cat( projCat, dNND)
catChild.selEventsFromID(dNND['aEqID_c'], repeats=True)
catParent.selEventsFromID(dNND['aEqID_p'], repeats=True)
print('before::: size of parent catalog', catParent.size(), 'size of offspring cat', catChild.size())

# ================================================================================
#                       bigger parent event pairs
# ================================================================================
# select event pairs with parent event larger than M_pt

sel = catParent.data['Mag'] >= dPar['M_pt']
catChild.selEventsFromID(dNND['aEqID_c'][sel], repeats=True)
catParent.selEventsFromID(dNND['aEqID_p'][sel], repeats=True)
print('after::: size of parent catalog', catParent.size(), 'size of offspring cat', catChild.size())

# ================================================================================
#                       calculate R/L
# ================================================================================
dsigma = 5e6
Mw = catParent.data['Mag']
M0 = 10 ** (1.5 * (Mw + 6.03))
l = 0.001 * (7 / 16 * M0 / dsigma) ** (1 / 3)
## hd
HD = dNND['aHD'][sel]
## Rl
aRl = HD / l

# ================================================================================
#                 select event paris in tails
# ================================================================================
sel_ev = aRl >= dPar['R_Lmin']
catChild.selEventsFromID(dNND['aEqID_c'][sel][sel_ev], repeats=True)
catParent.selEventsFromID(dNND['aEqID_p'][sel][sel_ev], repeats=True)
print('after::: size of parent catalog', catParent.size(), 'size of offspring cat', catChild.size())
aChild_lon = catChild.data['Lon']
aChild_lat = catChild.data['Lat']
aChild_mag = catChild.data['Mag']
aParent_lon = catParent.data['Lon']
aParent_lat = catParent.data['Lat']
aParent_mag = catParent.data['Mag']
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
        mymap.plot(x, y, marker='None', color='black', linewidth=1)
        x = []
        y = []

# labels = [left,right,top,bottom]
parallels = np.arange(dPar['minlat'], dPar['maxlat'], 2)
mymap.drawparallels(parallels, labels=[True, False, True, False])
meridians = np.arange(dPar['minlon'], dPar['maxlon'], 2)
mymap.drawmeridians(meridians, labels=[True, False, False, True])

# plot events
X_Child, Y_Child = mymap(aChild_lon, aChild_lat)
X_Parent, Y_Parent = mymap(aParent_lon,aParent_lat)
eqChildren = mymap.scatter(X_Child, Y_Child, aChild_mag, marker='o', c='blue', alpha=1)  # c=times,
eqParents = mymap.scatter(X_Parent, Y_Parent, aParent_mag, marker='o', c='red', alpha=1)  # c=times,

# plot parent-offspring pair
x = []
y = []
for x2,x1,y2,y1 in zip(X_Child,X_Parent,Y_Child,Y_Parent):
    x=[x1,x2]
    y=[y1,y2]
    if (x2-x1)**2+(y2-y1)**2 >= 10:
        print('***************',(x2-x1)**2+(y2-y1)**2)
    mymap.plot(x, y, marker='None', color='black', linewidth=1,alpha=0.2)
    x = []
    y = []

### plot landers and coso
catLC, ___,____,____=getLanders()
xC,yC=mymap(catLC.data['Lon'][1:],catLC.data['Lat'][1:])
xL,yL=mymap(catLC.data['Lon'][0],catLC.data['Lat'][0])
eqCoso = mymap.scatter(xC, yC, catLC.data['Mag'][1:], marker='s', c='orange', alpha=1)  # c=times,
eqLanders = mymap.scatter(xL, yL, catLC.data['Mag'][0], marker='*', c='red', alpha=1)  # c=times,


plt.show()
fig.savefig('plots/map_landers_coso_2',dpi=500)

"""
for curr_Mc,k in zip(dPar['aMc'],dPar['k']):
    # cut below current completeness
    eqCatMc.copy(eqCat)
    eqCatMc.selector(curr_Mc, None, 'Mag')
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
    print("haversine distance:", HD)
    print("rupture dimension", l)
    # ================================================================================
    #                       bigger parent event pairs
    # ================================================================================
    # select event pairs with parent event larger than M_pt
    '''
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
    '''
    # ================================================================================
    #                 select event paris in tails
    # ================================================================================
    sel_ev = aRl >= dPar['R_Lmin']
    catChild.selEventsFromID(dNND['aEqID_c'][sel_ev], repeats=True)
    catParent.selEventsFromID(dNND['aEqID_p'][sel_ev], repeats=True)
    print('after::: size of parent catalog', catParent.size(), 'size of offspring cat', catChild.size())
    aChild_lon = catChild.data['Lon']
    aChild_lat = catChild.data['Lat']
    aChile_mag = catChild.data['MAG']
    aParent_lon = catParent.data['Lon']
    aParent_lat = catParent.data['Lat']
    aParent_mag = catParent.data['MAG']
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
            mymap.plot(x, y, marker='None', color='black', linewidth=1)
            x = []
            y = []

    # labels = [left,right,top,bottom]
    parallels = np.arange(dPar['minlat'], dPar['maxlat'], 2)
    mymap.drawparallels(parallels, labels=[True, False, True, False])
    meridians = np.arange(dPar['minlon'], dPar['maxlon'], 2)
    mymap.drawmeridians(meridians, labels=[True, False, False, True])

    # plot events
    X_Child,Y_Child=mymap(aChild_lon,aChild_lat)
    eqs = mymap.scatter(X_Child, Y_Child, aChile_mag, marker='o',  vmin=0, vmax=maxtime, cmap=cmap, alpha=1) # c=times,
    plt.show()
"""