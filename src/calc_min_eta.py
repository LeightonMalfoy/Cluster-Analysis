
import numpy as np
import time
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import DistanceMetric

def haversine( lon1, lat1, lon2, lat2, **kwargs):
    """
    haversine formula implementation
    https://en.wikipedia.org/wiki/Great-circle_distance
    great circle distance between two points
    :input   lon1, lat1  - location of first set of points
             lon2, lat2  - loc. second set of points
                          - could be arrays or floating points


    :output  distance - great circle distance in kilometer
    """
    gR = 6378.137 # ~6370 # gR - Earth radius
    # convert to radians
    lon1 = lon1 * np.pi / 180
    lon2 = lon2 * np.pi / 180
    lat1 = lat1 * np.pi / 180
    lat2 = lat2 * np.pi / 180
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    distance = gR * c
    return distance

def calc_eta(eq1,eq2):
    """
    :param eq1: unit: decimal year
    :param eq2: unit: decimal year
    :return: distance defined as eta between two earthquakes ( 1 before 2)
    """
    dPar = {'D':1.6,
            'b':1.0,
            'm0':0, }
    time1, time2 = (eq1.data['Time'], eq2.data['Time'])
    lon1, lon2 = (eq1.data['lon'], eq2.data['lon'])
    lat1, lat2 = (eq1.data['lat'], eq2.data['lat'])
    m1 = eq1.data['Mag']
    hd = haversine(lon1,lat1,lon2,lat2)
    if time2>time1:
        eta = (time2-time1)*hd**dPar['D']*10**(dPar['b']*(m1-dPar['m0']))
    else:
        eta = 100000000
    return(eta)


def calc_min_eta(eqCat, dPar, algorithm='ball_tree'):
    """calculate min Nearest Neighbor Distance
    """
    t0 = time.time()
    aDist = np.zeros(eqCat.size()-1)
    time1,time2 = (eqCat.data['Time'],eqCat.data['Time'])
    lon1,lon2 = (eqCat.data['lon'],eqCat.data['lon'])
    lat1,lat2 = (eqCat.data['lat'],eqCat.data['lat'])
    m1 = eqCat.data['Mag']
    DistanceMetric.get_metric('pyfunc',kwargs=calc_eta)
    model = NearestNeighbors(n_neighbors=1, metric=calc_eta,algorithm=algorithm).fit(eqCat)

    # calculate distances knn method
    for i in np.arange(eqCat.size()):
        d_min_dist = model.kneighbors(eqCat)[0].min(axis=1)
        dist_arr[i, :] = d_min_dist
    print(time.time() - t0)
    return dist_arr