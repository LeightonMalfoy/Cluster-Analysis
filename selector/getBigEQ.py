"""
get the calalog containing one mainshock and its offsprings (1st generation aftershocks)

@author: Litong Huang, Peking Univ
"""

import sys
sys.path.append('/auto/home/lhuang/PycharmProjects/clustering-analysis-master')
from src.EqCat import *

def getMajor(curr_Mc,**kwargs):
    """
    :param curr_Mc: low-cut magnitude
    :return:
        - aCluster: list of clusters
          each cluster is a catalog contains a major shock and its offsprings
          'sel_p','sel_c' is used to select related NND data in aNND
        - aNND: nearest neighbor distance
          'aEqID_p': parent's ID of the pair
          'aEqID_c': child's ID of the pair
    """
    import numpy as np
    import os

    # ------------------------------my modules--------------------------------------
    import sys
    sys.path.append('/auto/home/lhuang/PycharmProjects/clustering-analysis-master')
    import src.data_utils as data_utils


    eqCat = EqCat()  # original catalog
    eqCatMc = EqCat()  # this catalog will be modified with each Mc iteration

    # =================================1==============================================
    #                            dir, file, params
    # ================================================================================
    data_dir = './data' # Todo: .. or .
    file_in = 'hs_1981_2018_all.mat'

    eqCat.loadMatBin(os.path.join(data_dir, file_in))
    eqCat.toCart_coordinates()
    eqCatMc.copy(eqCat)
    if 'min' in kwargs.keys():
        min = kwargs['min']
    else:
        min = 6.0
    if 'max' in kwargs.keys():
        max = kwargs['max']
    else:
        max = None
    eqCatMc.selectEvents(min,max,'Mag')

    # load nearest neighbor distances
    NND_file = './data/%s_NND_Mc_%.1f_HD.mat' % (file_in.split('.')[0], curr_Mc)# Todo: .. or .
    dNND = data_utils.loadmat(NND_file)  # ,  struct_as_record=True)

    aCluster = np.array([])
    for i in list(range(eqCatMc.size())):
        cat = EqCat()
        cat.copy(eqCat)
        sel_c = dNND['aEqID_p'] == eqCatMc.data['N'][i]
        sel_p = dNND['aEqID_c'] == eqCatMc.data['N'][i]
        sel = np.logical_or(sel_p,sel_c)
        cat.selEventsFromID(dNND['aEqID_c'][sel],repeats=True)
        cat.data['sel_p'] = sel_p
        cat.data['sel_c'] = sel_c
        aCluster=np.append(aCluster,cat)
        print("major earthquake:%.1f"%cat.data['Time'][0],"mag:",cat.data['Mag'][0])
    print("Total Ms: %d" % aCluster.shape[0])
    return aCluster,dNND

