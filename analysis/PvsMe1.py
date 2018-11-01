import numpy as np
import pylab as plt
import os, sys
import ipdb
from pylatex import Document, Package, Section, Figure, NoEscape, SubFigure
rootFolder = "" #"~/Documents/lab"
basefolder = rootFolder + "/homecentral/srao/Documents/code/mypybox"
print basefolder
sys.path.append(rootFolder + '/homecentral/srao/Documents/code/mypybox/utils')
from Print2Pdf import Print2Pdf
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import chdtrc as GammaQ
from scipy.signal import argrelextrema
import matplotlib.animation as animation
#import statsmodels.api as statsmodels
from matplotlib.ticker import FormatStrFormatter

def LoadFr(p, gamma, phi, mExt, mExtOne, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False):
    baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
    
    if nPop == 1:
    	baseFldr = baseFldr + 'onepop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if nPop == 2:
	if gamma >= .1 or gamma == 0:
	    baseFldr = baseFldr + 'twopop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
	else:
	    baseFldr = baseFldr + 'twopop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(T*1e-3), trNo)	    
    if IF_VERBOSE:
    	print baseFldr

    # if nPop == 1:
    # 	baseFldr = baseFldr + 'onepop/data/old/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    # if nPop == 2:
    # 	baseFldr = baseFldr + 'twopop/data/old/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    # if IF_VERBOSE:
#    print baseFldr
	
    filename = 'meanrates_theta%.6f_tr%s.txt'%(phi, trNo)
    print filename
    return np.loadtxt(baseFldr + filename)


def KappaVsM1AtPhi(kappaList, p=0, gamma=0, phi=0, mExt=0.075, mExtOne=0.075, rewireType='rand', N=10000, K=1000, nPop=2, T=1000, trNo=0, IF_PO_SORTED = False, sortedIdx = [], minRate = 0):
    m1E = np.empty((len(kappaList, )))
    validKappa = np.empty((len(kappaList, )))
    m1E[:] = np.nan
    validKappa[:] = np.nan
    mrMean = np.nan
    for kIdx, kappa in enumerate(kappaList):
	try:
            IF_VERBOSE = True
            if trNo == 0:
                print 'loading from fldr: ',
            #     IF_VERBOSE = True
	    mr = LoadFr(kappa, gamma, phi, mExt, mExtOne, trNo, T, N, K, nPop, IF_VERBOSE = IF_VERBOSE)
            mrMean = np.mean(mr[:N])
	    # if kappa == 16:
	    # 	ipdb.set_trace()

	    if not IF_PO_SORTED:
		m1E[kIdx] = M1Component(mr[:N])
	    else:
		mre = mr[:N]
		#m1E[kIdx] = M1Component(mre[sortedIdx[kIdx]])
                mask = mr[sortedIdx] > minRate
                m1E[kIdx] = M1Component(mr[sortedIdx[mask]])		
	    validKappa[kIdx] = kappa
	    # ipdb.set_trace()
	    print 'o',
	except IOError:
	    print 'x', 
	    #pass
	    #print 'kappa: ', kappa, ' no files!'
    sys.stdout.flush()	    
    return validKappa, m1E, mrMean



def KappaVsM1AtTr(kappaList, p=0, gamma=0, nPhis = 8, mExt=0.075, mExtOne=0.075, rewireType='rand', N=10000, K=1000, nPop=2, T=1000, trNo=0, IF_PO_SORTED = False, sortedIdx = [], minRate = 0):
    thetas = np.linspace(0, 180, nPhis, endpoint = False)
    m1E = np.zeros((nPhis, len(kappaList)))
    m0E = np.zeros((nPhis, len(kappaList)))    
    vldKappa = np.empty((nPhis, len(kappaList)))
    vldKappa[:] = np.nan
    # sortedIdx = []
    # if IF_PO_SORTED:
    # 	for kappa in kappaList:
    # 	    try:
    # 		po = GetPOofPop(p, gamma, mExt, 0.075, rewireType, 8, 0, N, K, nPop, T, IF_IN_RANGE = True, kappa = 1) # returns only for E neurons
    # 		sortedIdx.append(np.argsort(po))
    # 		# print sortedIdx[:10]
    # 	    except IOError:
    # 		pass
    for i, phi in enumerate(thetas):
	print 'phi =', phi
	vldKappa[i, :], m1E[i, :], m0E[i, :] = KappaVsM1AtPhi(kappaList, p, gamma, phi, mExt, mExtOne, rewireType, N, K, nPop, T, trNo, IF_PO_SORTED = IF_PO_SORTED, sortedIdx = sortedIdx, minRate = minRate)
    return np.nanmean(m1E, 0), np.nanmean(m0E, 0) #, vldKappa

def KappaVsM1(kappaList, nTrials = 10, p=0, gamma=0, nPhis = 8, mExt=0.075, mExtOne=0.075, rewireType='rand', N=10000, K=1000, nPop=2, T=1000, IF_PO_SORTED = False, sortedIdx = [], minRate = 0):
    m1E = np.zeros((nTrials + 2, len(kappaList)))
    m1E[0, :] = np.nan
    m0E = np.zeros((nTrials + 2, len(kappaList)))
    m0E[0, :] = np.nan    
    for trNo in range(0, nTrials + 2): #$ trNo 0 is always the CONTROL 
	print ''
	print '--' * 27
	print 'tr#: ', trNo
	print '--' * 27
        # m1E[trNo, :] = KappaVsM1AtTr(kappaList, p, gamma, nPhis, mExt, mExtOne, rewireType, N, K, nPop, T, trNo, IF_PO_SORTED, sortedIdx)
        m1E[trNo, :], m0E[trNo, :] = KappaVsM1AtTr(kappaList, p, gamma, nPhis, mExt, mExtOne, rewireType, N, K, nPop, T, trNo, IF_PO_SORTED, sortedIdx, minRate = minRate)
    m1EvsKappa = np.nanmean(m1E, 0)
    nValidTrials = np.sum(~np.isnan(m1E.mean(1)))
    if nValidTrials >= 2:
	m1EvsKappaSEM = np.nanstd(m1E, 0) / np.sqrt(float(nValidTrials))
    else:
	m1EvsKappaSEM = np.nan
    (_, caps, _) = plt.errorbar(kappaList, m1EvsKappa, fmt = 'o-', markersize = 3, yerr = m1EvsKappaSEM, lw = 0.8, elinewidth=0.8, label = r'$N = %s, K = %s$'%(N, K))
    for cap in caps:
        cap.set_markeredgewidth(0.8)
    plt.xlim(min(kappaList) - 1, max(kappaList) + 1)
    plt.xlabel(r'$\kappa$', fontsize = 20)
    plt.ylabel(r'$m_E^{(1)}$', fontsize = 20)
    plt.show()
    return m1EvsKappa, m1EvsKappaSEM
