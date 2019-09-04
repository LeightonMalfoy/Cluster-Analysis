"""
plot the histogram of a certain value's distribution
in a linear scale or in a logarithmic scale

"""

import matplotlib.pyplot as plt
import numpy as np

def plot_loglog(aBins,aHist,curr_Mc,name,**kwargs):
    plt.figure()
    ax=plt.subplot()
    ax.loglog(aBins, aHist, 'ko', mfc = 'none', mew = .2, label='Historical Rate')
    if 'aRate_bin' in kwargs.keys() and 'aRate' in kwargs.keys():
        ax.loglog(kwargs['aRate_bin'], kwargs['aRate'], 'r-', label='Earthquake rate,Mc=%.1f'%curr_Mc)
    if 'n' in kwargs.keys() and 'c' in kwargs.keys():
        sel = aBins>=0.5*kwargs['L']
        aFit = np.zeros(len(aBins[sel]))
        aFit = kwargs['c']*aBins[sel]**(kwargs['n'])
        ax.loglog(aBins[sel],aFit,'g-',lw=1.5,label='Decay of aftershock density')
    if 'L' in kwargs.keys():
        L = kwargs['L']
        #L = 10**(-3.22+0.69*kwargs['Mag'])/2
        ax.loglog([L, L], ax.get_ylim(), 'w-', lw=1.5)
        ax.plot([L, L], ax.get_ylim(), 'k--', lw=1.5, label='Rupture Radius: %.1fkm' % L)
    if 'q' in kwargs.keys() and 'd' in kwargs.keys():
        q = kwargs['q']
        d = kwargs['d']
        gamma = kwargs['gamma']
        mag = kwargs['Mag']
        #L = 10**(0.44*mag-2)
        L = 0.00199*10**(0.49*mag)
        sel = aBins>=10
        aFit = np.zeros(len(aBins[sel]))
        r=aBins[sel]
        R =1+(10/L)**(gamma+1)
        beta = q/d*(R**(d/(1+gamma)))/((R**(q/(1+gamma)))+q/d-1)
        aFit = beta*d*(r**gamma)/((L**(gamma+1))*((r/L)**(gamma+1)+1)**(1+d/(gamma+1)))
        ax.loglog(r,aFit,'b-',lw=1.5,label='Decay of aftershock density')
    #plt.loglog([3.1,3.1], [min(aRate),max(aRate)], 'b-', label='2.5')
    print(aBins[np.argmax(aHist)])
    ax.legend( loc='upper right',fontsize=8)
    ax.set_ylabel('Probability Density')
    #ax.set_xlim(None, np.log(40))
    plt.grid('on')
    # set the title with mainshock info
    title = name
    if 'YR' in kwargs.keys():
        title += ' YR:%d' % kwargs['YR']
    if 'Mag' in kwargs.keys():
        title += ' Mag:%.1f' % kwargs['Mag']
    plt.title(title)

    if 'plotkind' in kwargs.keys():
        kind = kwargs['plotkind']
        if kind[0:4] in ['lgHD','lgHD1','lgHD2','lgHD3','lgHD4']:
            xlabel = 'Haversine Distance (km)'
        elif kind == 'lgRL' or kind=='lgRL1':
            xlabel = 'R/L'
    else:
        kind = 'unknown'
    plt.xlabel('%s' % xlabel)
    plotFile = 'plots/h_%s_%s.png' % (kind, name)

    print('save plot', plotFile)
    plt.savefig(plotFile)
    plt.clf()
    return 0

def plot_cumlog(aBins,aHist,aCum,curr_Mc,name,**kwargs):
    plt.figure()
    ax=plt.subplot()
    ax.loglog(aBins, aCum, 'r-', label='CDF, Mc=%.1f' % curr_Mc)

    ax.set_ylabel('Probability Density')
    #ax.set_xlim(None, np.log(40))
    plt.grid('on')
    # set the title with mainshock info
    title = name
    sel = aBins >= 10
    aFit = np.zeros(len(aBins[sel]))
    L = 0.73118*10 ** (0.58 * kwargs['Mag']-3)
    r = aBins[sel]
    q = 0.35
    d = 1.2
    gamma = 0.6
    R = 1 + (10 / L) ** (gamma + 1)
    beta = q / d * (R ** (d / (1 + gamma))) / ((R ** (q / (1 + gamma))) + q / d - 1)
    aFit = beta * d * (r ** gamma) / ((L ** (gamma + 1)) * ((r / L) ** (gamma + 1) + 1) ** (1 + d / (gamma + 1)))
    anti_sel = aBins<10
    aFit = np.concatenate((aHist[anti_sel],aFit),axis=0)
    afitcum = np.cumsum(aFit)
    ax.loglog(r, afitcum[sel], 'b-', lw=1.5, label='Decay of aftershock density')
    if 'YR' in kwargs.keys():
        title += ' YR:%d' % kwargs['YR']
    if 'Mag' in kwargs.keys():
        title += ' Mag:%.1f' % kwargs['Mag']
    plt.title(title)
    if 'plotkind' in kwargs.keys():
        kind = kwargs['plotkind']
        if kind in ['lgHD','lgHD1','lgHD2','lgHD3','lgHD4','HDcum']:
            xlabel = 'Haversine Distance (km)'
        elif kind == 'lgRL' or kind=='lgRL1':
            xlabel = 'R/L'
    else:
        kind = 'unknown'
    plt.xlabel('%s' % xlabel)
    ax.legend( loc='lower right',fontsize=8)
    plotFile = 'plots/h_%s_%s.png' % (kind, name)

    print('save plot', plotFile)
    plt.savefig(plotFile)
    plt.clf()
    return 0


def plot_rate_pd(aBins,aHist,curr_Mc,name,binsize,**kwargs):
    """
    :param curr_Mc: magnitude completeness
    :param kwargs: keys of eqCat, aRate_bin,aRate,
    :return:
    """
    # plot the histogram
    plt.figure()
    ax = plt.subplot()
    ax.bar(aBins, aHist, width=0.8 * binsize, align='edge', color='.5', label='Mc = %.1f' % (curr_Mc))
    print(aBins[np.argmax(aHist)])
    ax.legend(loc='upper right')
    ax.set_ylabel('Probability Density')
    # plot the dividing line
    if 'threshold' in kwargs.keys():
        threshold = kwargs['threshold']
        ax.plot([threshold, threshold], ax.get_ylim(), 'w-', lw=2)
        ax.plot([threshold, threshold], ax.get_ylim(), 'k--', lw=2)
    if 'L' in kwargs.keys():
        L = kwargs['L']
        ax.plot([L, L], ax.get_ylim(), 'w-', lw=1.5)
        ax.plot([L, L], ax.get_ylim(), 'k--', lw=1.5,label='Rupture Extent: %.1fkm'%L)
    # plot the linear rate
    if 'aRate_bin' in kwargs.keys() and 'aRate' in kwargs.keys():
        ax.plot(kwargs['aRate_bin'],kwargs['aRate'],'r-',lw=1,label='Linear rate')
        ax.legend(loc='upper right')
    # set the title with mainshock info
    title = name
    if 'YR' in kwargs.keys():
        title += ' YR:%d' % kwargs['YR']
    if 'Mag' in kwargs.keys():
        title += ' Mag:%.1f' % kwargs['Mag']
    ax.set_title(title)
    ax.grid('on')

    #ax.set_xlim(0, 200)
    plt.show()
    if 'plotkind' in kwargs.keys():
        kind = kwargs['plotkind']
    else:
        kind = 'unknown'
    ax.set_xlabel('%s' % kind)
    plotFile = 'plots/h_%s_%s.png' % (kind, name)

    print('save plot', plotFile)
    plt.savefig(plotFile)
    plt.clf()
    return 0

