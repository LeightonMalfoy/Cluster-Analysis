"""
plot major earthquakes and their offsprings':
    1. real distance distribution (histogram and linear rate)
    2. lg R -lg T distribution (bg: southern california distribution)
    3. NND histogram
    4. R/L distribution
"""
import sys
sys.path.append('/auto/home/lhuang/PycharmProjects/clustering-analysis-master')
import matplotlib.pyplot as plt
import numpy as np
from src.data_utils import density_2D

def plot_HD(filename,aHD,majorEQ,binsize,Mc):
    aBins = np.arange(min(aHD),max(aHD),binsize)
    aHist, aBins = np.histogram(aHD,aBins)
    aHist, aBins = np.array(list(zip(*zip(aHist, aBins))))  # cut to same length
    # quantity -> density
    aHist /= binsize
    # to pdf (prob. density)
    aHist /= len(aHD)
    # plot
    fig,ax = plt.subplots()
    plt.title("Year:%d Mag:%.1f"% (majorEQ.data['Time'],majorEQ.data['Mag']))
    ax.bar( aBins, aHist, width =.8*binsize, align = 'edge', color = '.5', label = 'Mc = %.1f'%( Mc))
    ax.legend(loc='upper right')
    ax.set_xlabel('Haversine Distance')
    ax.set_ylabel('Historical Density')
    ax.grid('on')
    ax.set_xlim(0,500)
    plt.savefig(filename,dpi=500)
    plt.show()
    plt.clf()

def plot_HD1(filename,aHD,aHD1,majorEQ,binsize,Mc):
    aBins = np.arange(min(aHD),max(aHD),binsize)
    aBins1 = np.arange(min(aHD), max(aHD), binsize)
    aHist, aBins = np.histogram(aHD,aBins)
    aHist1, aBins1 = np.histogram(aHD1,aBins1)
    aHist, aBins = np.array(list(zip(*zip(aHist, aBins))))  # cut to same length
    aHist1, aBins1 = np.array(list(zip(*zip(aHist1, aBins1))))  # cut to same length
    # quantity -> density
    aHist /= binsize
    aHist1 /= binsize
    # to pdf (prob. density)
    aHist /= len(aHD)
    aHist1 /= len(aHD)
    # plot
    fig,ax = plt.subplots()
    plt.title("Year:%d Mag:%.1f MC:%.1f"% (majorEQ.data['Time'],majorEQ.data['Mag'],Mc))
    ax.bar( aBins, aHist, width =.8*binsize, align = 'edge', color = '.5', label = 'All Event Pairs')
    ax.bar( aBins1, aHist1, width =.8*binsize, align = 'edge', color = 'r', label = 'Clustered Event Pairs')
    ax.legend(loc='upper right')
    ax.set_xlabel('Haversine Distance')
    ax.set_ylabel('Historical Density')
    ax.grid('on')
    ax.set_xlim(0,500)
    plt.savefig(filename,dpi=500)
    plt.show()
    plt.clf()

def plot_eta(filename,aEta,majorEQ,binsize,Mc):
    aBins = np.arange(min(aEta),max(aEta),binsize)
    aHist, aBins = np.histogram(aEta,aBins)
    aHist, aBins = np.array(list(zip(*zip(aHist, aBins))))  # cut to same length
    # quantity -> density
    aHist /= binsize
    # to pdf (prob. density)
    aHist /= len(aEta)
    # plot
    fig,ax = plt.subplots()
    plt.title("Year:%d Mag:%.1f"% (majorEQ.data['Time'],majorEQ.data['Mag']))
    ax.bar( aBins, aHist, width =.8*binsize, align = 'edge', color = '.5', label = 'Mc = %.1f'%( Mc))
    ax.legend(loc='upper right')
    ax.set_xlabel('lg(NND)')
    ax.set_ylabel('Historical Density')
    ax.grid('on')
    #ax.set_xlim(0,500)
    plt.savefig(filename,dpi=500)
    plt.show()
    plt.clf()

def plot_r_tau_density(filename,aTau,aR,dPar,majorEQ):
    a_Tbin = np.arange(dPar['Tmin'], dPar['Tmax'] + 2 * dPar['binx'], dPar['binx'])
    a_Rbin = np.arange(dPar['Rmin'], dPar['Rmax'] + 2 * dPar['biny'], dPar['biny'])
    XX, YY, ZZ = density_2D(aTau, aR, a_Tbin, a_Rbin, sigma=dPar['sigma'])
    plt.figure(1, figsize=(8, 10))
    ax = plt.subplot(111)
    normZZ = ZZ * (dPar['binx'] * dPar['biny'] * len(aTau))
    plot1 = ax.pcolormesh(XX, YY, normZZ, cmap=dPar['cmap'])
    cbar = plt.colorbar(plot1, orientation='horizontal', shrink=.5, aspect=20, )
    # plot eta_0 to divide clustered and background mode
    ax.plot([dPar['Tmin'], dPar['Tmax']], -np.array([dPar['Tmin'], dPar['Tmax']]) + dPar['eta_0'], '-', lw=1.5,
            color='w')
    ax.plot([dPar['Tmin'], dPar['Tmax']], -np.array([dPar['Tmin'], dPar['Tmax']]) + dPar['eta_0'], '--', lw=1.5,
            color='.5')
    cbar.set_label( 'Number of Event Pairs',labelpad=-40)
    ax.set_xlabel( 'Rescaled Time')
    ax.set_ylabel( 'Rescaled Distance')
    ax.set_title("Year:%d Mag:%.1f/nNearest Neighbor Pairs in R-T"% (majorEQ.data['Time'],majorEQ.data['Mag']))
    plt.savefig(filename,dpi=500)
    plt.clf()