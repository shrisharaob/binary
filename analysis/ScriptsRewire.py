''' USAGE: [dbName, NetworkType, K, NE, NI, thetaSig] = DefaultArgs(sys.argv[1:], ['', 'ori', 1000, 10000, 10000, 0.5]) '''

basefolder = "/homecentral/srao/Documents/code/mypybox"
import numpy as np
import code, sys, os
import ipdb
import pylab as plt
sys.path.append(basefolder)
import Keyboard as kb
from multiprocessing import Pool
from functools import partial 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
sys.path.append(basefolder + "/nda/spkStats")
sys.path.append(basefolder + "/utils")
from DefaultArgs import DefaultArgs
from reportfig import ReportFig
from Print2Pdf import Print2Pdf
import GetPO


rootFolder = ''

def GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2):
    if rewireType == 'rand':
	tag = ''
    if rewireType == 'exp':
	tag = '1'

    rootFolder = ''
    baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
    if nPop == 1:
    	baseFldr = baseFldr + 'onepop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if nPop == 2:
	if gamma >= .1 or gamma == 0:
	    baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
	else:
	    baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/p%sgamma/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(T*1e-3), trNo)	        
    return baseFldr

def LoadFr(p, gamma, phi, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False):
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop)    
    if IF_VERBOSE:
    	print baseFldr
    filename = 'meanrates_theta%.6f_tr%s.txt'%(phi, trNo)
    print filename
    return np.loadtxt(baseFldr + filename)

def GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000):
    NE = N
    NI = N
    tc = np.zeros((NE + NI, nPhis))
    tc[:] = np.nan
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    for i, iPhi in enumerate(phis):
	print i, iPhi
	try:
	    if i == 0:
		print 'loading from fldr: ',	    
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = True)
	    else:
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = False)
	    if(len(fr) == 1):
		if(np.isnan(fr)):
		    print 'file not found!'
	    tc[:, i] = fr
	except IOError:
	    print 'file not found!'
	    raise SystemExit
    return tc

def FollowNeurons(neuronsList, p, gamma, nPhis, mExt, mExtOne, rewireType, trList = [0], N = 10000, K = 1000, nPop = 2, T = 1000, IF_LEGEND = True):
        tc = [[] for x in xrange(len(trList))]
	preferredOri = [[] for x in xrange(len(trList))]
	theta = np.linspace(0, 180, nPhis, endpoint = False)
	for kk, trNo in enumerate(trList):
	    tmp = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T)
	    tc[kk].append(tmp)
	    preferredOri[kk].append(POofPopulation(tc[kk][0][:N], IF_IN_RANGE = True))
        # ipdb.set_trace()
	plt.ion()
	for neuronIdx in neuronsList:
	    for kk, trNo in enumerate(trList):
		osi = OSI(tc[kk][0][neuronIdx, :], theta)
		plt.plot(theta, tc[kk][0][neuronIdx, :], 's-', label = 'STP:%s OSI:%.4s PO:%.5s'%(trNo, osi, preferredOri[kk][0][neuronIdx]))
		plt.title('neuron#%s'%(neuronIdx))
		plt.xlabel('Stimulus (deg)')
		plt.ylabel(r'$m(\phi)$')
	    if IF_LEGEND:
		plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 8})
	    plt.waitforbuttonpress()
	    plt.clf()

def FollowNeuron(neuronIdx, p, gamma, nPhis, mExt, mExtOne, rewireType, trList = [0], N = 10000, K = 1000, nPop = 2, T = 1000, IF_LEGEND = True):
        tc = [[] for x in xrange(len(trList))]
	preferedOri = [[] for x in xrange(N)]
	theta = np.linspace(0, 180, nPhis, endpoint = False)
	for trNo in trList:
	    tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T)
	    preferredOri = POofPopulation(tc[:N], IF_IN_RANGE = True)
	    osi = OSI(tc[neuronIdx, :], theta)
	    plt.plot(theta, tc[neuronIdx, :], 's-', label = 'STP:%s OSI:%.4s PO:%.5s'%(trNo, osi, preferredOri[neuronIdx]))
	plt.title('neuron#%s'%(neuronIdx))
	plt.xlabel('Stimulus (deg)')
	plt.ylabel(r'$m(\phi)$')
        if IF_LEGEND:
	    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 8})


def GetPhase(firingRate, atTheta, IF_IN_RANGE = False):
    out = np.nan
    zk = np.dot(firingRate, np.exp(2j * atTheta * np.pi / 180))
    out = np.angle(zk) * 180.0 / np.pi
    if IF_IN_RANGE:
	if(out < 0):
	    out += 360
    return out * 0.5

def POofPopulation(tc, theta = np.arange(0.0, 180.0, 22.5), IF_IN_RANGE = False):
    # return value in degrees
    nNeurons, nAngles = tc.shape
    theta = np.linspace(0, 180, nAngles, endpoint = False)
    po = np.zeros((nNeurons, ))
    for kNeuron in np.arange(nNeurons):
        po[kNeuron] = GetPhase(tc[kNeuron, :], theta, IF_IN_RANGE)
    return po 

def GetPOofPop(p, gamma, mExt, mExtOne, nPhis = 8, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_IN_RANGE = True):
    nNeurons = N
    thetas= np.linspace(0, np.pi, nNeurons, endpoint = False)
    tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, trNo, N, K, nPop, T)
    prefferedOri = POofPopulation(tc[:N], IF_IN_RANGE = True) * np.pi / 180.0
    return prefferedOri

def OSI(firingRate, atTheta):
    out = np.nan
    zk = np.dot(firingRate, np.exp(2j * atTheta * np.pi / 180))
    # denum = np.absolute(np.dot(firingRate, np.exp(2j * np.zeros(len(atTheta)))))
    if(firingRate.mean() > 0.0):
        out = np.absolute(zk) / np.sum(firingRate)
    return out

def OSIOfPop(firingRates, atThetas):
    # thetas in degrees
    nNeurons, nThetas = firingRates.shape
    out = np.zeros((nNeurons, ))
    for i in range(nNeurons):
        out[i] = OSI(firingRates[i , :], atThetas)
    return out

def POofPopulation(tc, theta = np.arange(0.0, 180.0, 22.5), IF_IN_RANGE = False):
    # return value in degrees
    nNeurons, nAngles = tc.shape
    theta = np.linspace(0, 180, nAngles, endpoint = False)
    po = np.zeros((nNeurons, ))
    for kNeuron in np.arange(nNeurons):
        po[kNeuron] = GetPhase(tc[kNeuron, :], theta, IF_IN_RANGE)
    return po 

def PltOSIHist(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_NEW_FIG = True, color = 'k'):
    NE = N
    NI = N
    tc = np.zeros((NE + NI, nPhis))
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    if IF_NEW_FIG:
	plt.figure()
    for i, iPhi in enumerate(phis):
	print i, iPhi
	try:
	    if i == 0:
		print 'loading from fldr: ',	    
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = True)
	    else:
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = False)
	    if(len(fr) == 1):
		if(np.isnan(fr)):
		    print 'file not found!'
	    tc[:, i] = fr
	except IOError:
	    print 'file not found!'
    osi = OSIOfPop(tc[:NE, :], phis)
    print "K = ", K, ", osi simulation: ", np.nanmean(osi)
    # plt.xlabel(r"$\mathrm{OSI} \,\,\,\,  (m_{E, i}^{(1)})$")
    plt.xlabel('OSI', fontsize = 12)    
    plt.ylabel('Density', fontsize = 12)
    plt.hist(osi[~np.isnan(osi)], 27, normed = 1, histtype = 'step', label = r'$STP: %s$'%(trNo, ), color = color, lw = 1)    
    # plt.hist(osi[~np.isnan(osi)], 27, normed = 1, histtype = 'step', label = r'$p = %s,\,\gamma = %s$'%(p, gamma, ), color = color, lw = 1)
    plt.xlim(0, 1)
    plt.gca().set_xticks([0, 0.5, 1])
    _, ymax = plt.ylim()
    plt.gca().set_yticks([0, np.ceil(ymax)])    
    # plt.title(r'$N = %s,\, K = %s,\, m_0^{(0)} = %s,\, m_0^{(1)} = %s$'%(NE, K, mExt, mExtOne))
    plt.title(r'$m_0^{(0)} = %s, \,m_0^{(1)} = %s $'%(mExt, mExtOne))
    _, ymax = plt.ylim()
    plt.vlines(np.nanmean(osi), 0, ymax, lw = 1, color = color)
    print "mean OSI = ", np.nanmean(osi)
    # print "MEAN OSI = ", np.nanmean(osi)
    return osi

def CompareOSIHist(pList, gList, nPhis, mExt, mExtOneList, trList, rewireType, N = 10000, KList = [1000], nPop = 2, T = 1000, IF_NEW_FIG = True, clrCntr = 0, filename = '', IF_LEGEND = True, legendTxt = ''):
    if IF_NEW_FIG:
	plt.figure()
    colors = [plt.cm.Dark2(i) for i in np.linspace(0, 1, 1 + clrCntr + len(pList) * len(gList) * len(mExtOneList) * len(trList), endpoint = False)]
    meanOSI = []
    for mExtOne in mExtOneList:
	for trNo in trList:
	    for p in pList:
		for gamma in gList:
		    for K in KList:
			try:
			    print trNo
			    tmposi = PltOSIHist(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo = trNo, IF_NEW_FIG = False, color = colors[clrCntr], T=T, K=K)
			    meanOSI.append(np.nanmean(tmposi))
			    clrCntr += 1
			except IOError:
			    print "p = ", p, " gamma = ", gamma, " trial# ", trNo, " file not found"
    # plt.gca().legend(bbox_to_anchor = (1.1, 1.5))
    if IF_LEGEND:
	plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 8})

    plt.gca().set_position([0.15, 0.15, .65, .65])
    if nPop == 2:
	# plt.savefig("./figs/twopop/compareOSI_.png")
	# plt.savefig("./figs/twopop/compareOSI_"+filename + '.png')
	paperSize = [4, 3]
	# ipdb.set_trace()
	filename = filename + "p%sg%s"%(p, gamma)
        Print2Pdf(plt.gcf(),  "./figs/twopop/compareOSI_"+filename,  paperSize, figFormat='png', labelFontsize = 10, tickFontsize=8, titleSize = 10.0, IF_ADJUST_POSITION = True, axPosition = [0.14, 0.14, .7, .7])
    print '--'*26
    print 'pc change in mean OSI = ', 100 * (meanOSI[-1] - meanOSI[0]) / meanOSI[0]
    print '--'*26
    
    

    
