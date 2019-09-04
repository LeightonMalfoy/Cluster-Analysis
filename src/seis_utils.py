#python2.7
"""

 utility functions to analyze earthquake data sets

"""
from __future__ import division
import numpy as np

#===================================================================================
#                         rate calculation
#===================================================================================
def eqRate( aValue, k_win,**kwargs):
    """
    smoothed rate from overlapping sample windows normalized by delta_value
    :param aValue:
    :param k_win: window size
    :param kwargs: control k decrease
    :return:
    """

    aValue = np.sort(aValue)
    aS          = np.arange( int(0.5*k_win), int(aValue.shape[0]-0.5*k_win), 1)
    aBin, aRate = np.array([]),np.array([])
    iS = 0
    for s in aS:
        i1, i2 = int(s-0.5*k_win), int(s+0.5*k_win)
        if i2 >= aValue.shape[0]:
            break
        aBin = np.append(aBin,0.5*( aValue[i1]+aValue[i2]))
        aRate = np.append(aRate, k_win/( aValue[i2]-aValue[i1]))
        if 'step' in kwargs.keys():
            step = -1*kwargs['step']
        else:
            step = -6
        if 'minK' in kwargs.keys() and 'dividing' in kwargs.keys():
            if (aBin[iS] >= kwargs['dividing']) & (k_win >= kwargs['minK']):
                k_win -= step
        iS += 1
    return aBin, aRate

def eqlgRate( algValue, k_bin,**kwargs):
    # calculate the earthquake rate in a logrithmic scale
    lgBinsize = (np.amax(algValue)-np.amin(algValue))/k_bin
    algBins = np.arange(np.amin(algValue),np.amax(algValue),lgBinsize)
    algHist,algBins = np.histogram(algValue,algBins)
    deltaV=10**algBins[1:]-10**algBins[0:len(algBins)-1]
    algRate = algHist/deltaV
    algRate /= len(algValue)
    algRate,algBins = np.array(list(zip(*zip(algRate,algBins))))
    aBins = 10**algBins
    return aBins,algRate

#===================================================================================
#                         time and distance
#===================================================================================
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

def dateTime2decYr( datetime_in, **kwargs ):
    """
    input: datetime_in = array containing time columns year - second
                   out = date in decimal year
                   
    """
    import datetime
    import calendar
    try:
        o_dt = datetime.datetime( int( datetime_in[0] ), int( datetime_in[1] ), int( datetime_in[2] ), int( datetime_in[3] ), int( datetime_in[4] ), int(datetime_in[5]))
    except:
        error_msg = "datetime array not valid - %s; check if date and time is correct, e.g. no SC > 60.." % datetime_in
        raise ValueError
        raise error_msg
    time_sc = o_dt.hour*3600 + o_dt.minute*60 + o_dt.second    
    # get no. of day within current year between 0 to 364 and ad time in seconds
    dayOfYear_seconds = ( o_dt.timetuple().tm_yday - 1 ) * 86400.0 + time_sc
    if calendar.isleap( o_dt.year):
        year_fraction = dayOfYear_seconds / ( 86400.0 * 366 )
    else:
        year_fraction = dayOfYear_seconds / ( 86400.0 * 365 )
    # dec year = current year + day_time (in dec year)
    return o_dt.year + year_fraction

# def dateTime2decYr_old( YR, MO, DY, HR, MN, SC):
#     """
#     - convert date time to decimal year
#     :param YR: - int or arrays
#     :param MO:
#     :param DY:
#     :param HR:
#     :param MN:
#     :param SC:
#     :return:
#     """
#     nDays = 365.25
#     return YR + (MO-1)/12 + (DY-1)/nDays  + HR/(nDays*24) + MN/(nDays*24*60) + SC/(nDays*24*3600)

def area_poly( aX, aY):
    """
    use:

    A = 0.1*abs( (x1*y2 + x2*y3 + xn-1*yn + xn*y1) - (y1*x2 + y2*x3 + ... + yn-1*xn + yn*x1))
    :param aX: - x-coordinates of all vertices
    :param aY: - y-coordinates of all vertices
    :return: A - area of polygon
    """
    #sumVert1 = (aX[0:-1]*aY[1::]).sum()+aX[-1]*aY[0]
    # or:
    sumVert1  = np.dot( aX[0:-1], aY[1::])+aX[-1]*aY[0]
    #sumVert2 = (aY[0:-1]*aX[1::]).sum()+aY[-1]*aX[0]
    # or:
    sumVert2  = np.dot(aY[0:-1], aX[1::])+aY[-1]*aX[0]
    #sum = (aX[0:-1]*aY[1::] - aY[0:-1]*aX[1::]).sum() + (aX[-1]*aY[0]-aY[-1]*aX[0])
    return 0.5*abs( sumVert1 - sumVert2)