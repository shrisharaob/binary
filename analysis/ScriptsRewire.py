''' USAGE: [dbName, NetworkType, K, NE, NI, thetaSig] = DefaultArgs(sys.argv[1:], ['', 'ori', 1000, 10000, 10000, 0.5]) '''

basefolder = "/homecentral/srao/Documents/code/mypybox"
import numpy as np
import code, sys, os
import ipdb
import pylab as plt
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
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
from scipy.optimize import curve_fit
sys.path.append('/homecentral/srao/Documents/code/binary/c/twopop')
import RewireNetwork as rw
rootFolder = ''
#from pylatex import Document, Package, Section, Figure, NoEscape, SubFigure

def GetBaseFolderOld(p, gamma, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, kappa = 1):
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
        

    print baseFldr	    
    
    return baseFldr

def GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, kappa = 1):
    # ipdb.set_trace()
    tag = ''
    if kappa == -1:
	baseFldr = GetBaseFolderOld(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    else:
	if rewireType == 'rand':
	    tag = ''
	if rewireType == 'exp':
	    tag = '1'
	if rewireType == 'decay':
	    tag = '2'
	if rewireType == 'twosteps':
	    tag = '3'	    
	rootFolder = ''
	baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
	if nPop == 1:
	    baseFldr = baseFldr + 'onepop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
	# ipdb.set_trace()	
	if nPop == 2:
	    if T < 1000:
		baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), trNo)
	    else:
		if gamma >= .1 or gamma == 0:
		    baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
		else:
		    baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(T*1e-3), trNo)
    if rewireType == 'cntrl':
        baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
        baseFldr = baseFldr + 'twopop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    print baseFldr
    return baseFldr

def LoadFr(p, gamma, phi, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False, kappa = 0):
    #ipdb.set_trace()
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    filename = 'meanrates_theta%.6f_tr%s.txt'%(phi, trNo)
    if IF_VERBOSE:
    	print baseFldr
	print filename
    return np.loadtxt(baseFldr + filename)

def LoadFrChnk(chnk, phi, rewireType, p=0, gamma = 0, mExt = 0.075, mExtOne = 0.075, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False, kappa = 0):
    # ipdb.set_trace()
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    filename = 'meanrates_theta%.6f_tr%s_chnk%s.txt'%(phi, trNo, chnk)
    if IF_VERBOSE:
    	print baseFldr
	print filename
    return np.loadtxt(baseFldr + filename)
    

def LoadFFInput(p, gamma, phi, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False, kappa = 1):
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    if IF_VERBOSE:
    	print baseFldr
    filename = 'meanFFinput_theta%.6f_tr%s_last.txt'%(phi, trNo)
    print filename
    return np.loadtxt(baseFldr + filename)

def LoadEEInput(p, gamma, phi, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False, kappa = 1):
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    if IF_VERBOSE:
    	print baseFldr
    filename = 'meaninput_E_theta%.6f_tr%s_last.txt'%(phi, trNo)
    filenameOld = 'meaninput_EE_theta%.6f_tr%s_last.txt'%(phi, trNo)    
    print filename
    try:
	return np.loadtxt(baseFldr + filename)
    except IOError:
	return np.loadtxt(baseFldr + filenameOld)

def LoadEIInput(p, gamma, phi, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False, kappa = 1):
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    if IF_VERBOSE:
    	print baseFldr
    filename = 'meaninput_I_theta%.6f_tr%s_last.txt'%(phi, trNo)
    filenameOld = 'meaninput_EI_theta%.6f_tr%s_last.txt'%(phi, trNo)        
    print filename
    try:
	return np.loadtxt(baseFldr + filename)
    except IOError:
	return np.loadtxt(baseFldr + filenameOld)

def LoadNetInput(p, gamma, phi, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False, kappa = 1):
    # ipdb.set_trace()
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    if IF_VERBOSE:
    	print baseFldr
    filename = 'meaninput_theta%.6f_tr%s_last.txt'%(phi, trNo)
    print filename
    # ipdb.set_trace()
    return np.loadtxt(baseFldr + filename)

def CircularCorrCoeff(x, y):
    # x and y must be in radians
    n = x.size
    nX = x.size
    nY = y.size
    if nX != nY:
    	n = np.max([nX, nY])
    if(nX != n):
    	x = np.ones((n, )) * x
    if(nY != n):
    	y = np.ones((n, )) * y

    numerator = 0
    for i in range(n - 1):
	for j in range(i + 1, n):
	    numerator += np.sin(2.0 * (x[i] - x[j])) * np.sin(2.0 * (y[i] - y[j]))
    denom1 = 0
    denom2 = 0
    for i in range(n - 1):
	for j in range(i + 1, n):
	    denom1 += np.sin(2.0 * (x[i] - x[j]))**2
	    denom2 += np.sin(2.0 * (y[i] - y[j]))**2
    denom = np.sqrt(denom1 * denom2)
    return numerator / denom
    
def GetInputTuningCurves(p, gamma = 0, nPhis = 8, mExt = .075, mExtOne = .075, rewireType = 'rand', trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 1, inputType = 'net'):
    NE = N
    NI = N

    #tc = np.zeros((NE, nPhis))    
    # 
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    if inputType == 'net':
	LoadFunc = LoadNetInput
    elif inputType == 'FF':
	LoadFunc = LoadFFInput
    elif inputType == 'EE':
        # tc = np.zeros((NE, nPhis))	
	LoadFunc = LoadEEInput
    elif inputType == 'EI':
        # tc = np.zeros((NE, nPhis))		
	LoadFunc = LoadEIInput
    tc = np.zeros((nPop*N, nPhis))
    tc[:] = np.nan
    for i, iPhi in enumerate(phis):
	print i, iPhi
	try:
	    if i == 0:
		print 'loading from fldr: ',	    
		fr = LoadFunc(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = True, kappa = kappa)
	    else:
		fr = LoadFunc(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = False, kappa = kappa)
	    if(len(fr) == 1):
		if(np.isnan(fr)):
		    print 'file not found!'
	    if i == 0:
		tc = np.zeros((fr.size, nPhis))
		tc[:] = np.nan
	    tc[:, i] = fr
	except IOError:
	    print 'file not found!'
	    # raise SystemExit
    # ipdb.set_trace()
    return tc

def GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 0, IF_SUCCESS = False):
    NE = N
    NI = N
    tc = np.zeros((NE + NI, nPhis))
    tc[:] = np.nan
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    IF_FILE_LOADED = False
    for i, iPhi in enumerate(phis):
	print i, iPhi
	try:
	    if i == 0:
		print 'loading from fldr: ',	    
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = True, kappa = kappa)
	    else:
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = False, kappa = kappa)
	    if(len(fr) == 1):
		if(np.isnan(fr)):
		    print 'file not found!'
	    tc[:, i] = fr
            IF_FILE_LOADED = True
	except IOError:
	    print 'file not found!'
	    # raise SystemExit
    if IF_SUCCESS:
        return tc, IF_FILE_LOADED
    else:
        return tc

def GetTuningCurvesChnk(chnk, rewireType, p = 0, gamma = 0, nPhis = 8, mExt = 0.075, mExtOne = 0.075, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 0, IF_SUCCESS = False):
    NE = N
    NI = N
    tc = np.zeros((NE + NI, nPhis))
    tc[:] = np.nan
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    IF_FILE_LOADED = False
    for i, iPhi in enumerate(phis):
	print i, iPhi
	try:
	    if i == 0:
		print 'loading from fldr: ',
		fr = LoadFrChnk(chnk, iPhi, rewireType, p, gamma, mExt, mExtOne, trNo, T, NE, K, nPop, IF_VERBOSE = True, kappa = kappa)
	    else:
		fr = LoadFrChnk(chnk, iPhi, rewireType, p, gamma, mExt, mExtOne, trNo, T, NE, K, nPop, IF_VERBOSE = False, kappa = kappa)
	    if(len(fr) == 1):
		if(np.isnan(fr)):
		    print 'file not found!'
	    tc[:, i] = fr
            IF_FILE_LOADED = True
	except IOError:
	    print 'file not found!'
	    # raise SystemExit
    if IF_SUCCESS:
        return tc, IF_FILE_LOADED
    else:
        return tc


        
def PlotInOutTuningCurve(nNeurons, p = 0, gamma = 0, nPhis = 8, mExt = .075, mExtOne = .075, rewireType = 'rand', trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 1):
    tcIn = GetInputTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
    tcOut = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
    neuronIdx = np.random.choice(xrange(N), nNeurons)
    theta = np.linspace(0, 180, nPhis, endpoint = False)
    poIn = POofPopulation(tcIn, theta, IF_IN_RANGE = True) * np.pi / 180
    poOut = POofPopulation(tcOut, theta, IF_IN_RANGE = True) * np.pi / 180
    # plt.figure()
    # plt.plot(poIn[:N], poOut[:N], 'k.')
    # plt.show()
    # ipdb.set_trace()
    # plt.plot(np.cos(2.0 * (poIn[:N] - poOut[:N])))
    print '-' * 10
    print 'avg = ', CircularCorrCoeff(poIn[:N], poOut[:N])
    print '-' * 10
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    plt.ion()
    
    for i in neuronIdx:
	ax1.plot(theta, tcIn[i, :], 'g.-', label = '')
	ax2.plot(theta, tcOut[i, :], 'ko-', markerfacecolor = 'None', label = r'')
	ax1.set_ylabel(r'$u(\phi)$', fontsize = 16)
	ax1.tick_params('y', colors='g')
	ax2.set_ylabel(r'$m(\phi)$', fontsize = 16)
	ymin1, ymax1 = ax1.get_ylim()
        ymin2, ymax2 = ax2.get_ylim()
	ymin = np.min([ymin1, ymin2])
	ymax = np.max([ymax1, ymax2])
	ax1.vlines(poIn[i], ymin1, ymax1, color = 'g')
	ax2.vlines(poOut[i], ymin2, ymax2, color = 'k')	
        # plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 10})
	ax1.set_title('neuron#%s'%(i))
	ax1.set_xlabel(r'$\phi$')
	fig.waitforbuttonpress()
	ax1.clear()
	ax2.clear()

def PlotInOutPOCorr(nRewireSteps, kappa = 1, p = 0, gamma = 0, nPhis = 8, mExt = .075, mExtOne = .075, rewireType = 'rand', trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_COMPUTE = True):
    outE = []; outI = []
    validTr = []
    if IF_COMPUTE:
	for trNo in range(nRewireSteps):
	    try:
		tcIn = GetInputTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
		tcOut = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
		theta = np.linspace(0, 180, nPhis, endpoint = False)
		poIn = POofPopulation(tcIn, theta, IF_IN_RANGE = True)
		poOut = POofPopulation(tcOut, theta, IF_IN_RANGE = True)
		print 'computing ccc... ',
		sys.stdout.flush()
		outE.append(CircularCorrCoeff(poIn[:N], poOut[:N])) #np.nanmean(np.cos(2.0 * (poIn[:N] - poOut[:N]))))
		print 'done'
		# outI.append(CircularCorrCoeff(poIn[:N], poOut[:N])) #np.nanmean(np.cos(2.0 * (poIn[N:] - poOut[N:]))))
		validTr.append(trNo)
	    except IOError:
		print ''
	ipdb.set_trace()
	# np.save('InOutCorr_p%sg%sk%s'%(int(100 * p), int(10 * gamma), int(10 * kappa)), [validTr, outE]))
	filename = './data/twopop/InOutCorr_p%sg%sk%s'%(int(100 * p), int(10 * gamma), int(10 * kappa))
        np.save(filename, [validTr, outE])
	# plt.figure()
	# plt.plot(validTr, outE, 'ko-')
	# plt.xlabel('rewire step')
	# plt.ylabel(r'$\langle \cos 2 [ \phi_{j}^{po} - \theta_{j, IN}^{po} ]   \rangle_j$')
	# plt.plot(outI, 'ro-')
	filepath = './figs/twopop/InOutCorr_p%sg%s'%(int(100 * p), int(10 * gamma))
	ProcessFigure(plt.gcf(), filepath, True)

# sr.ComputeInOutPOCorrParallelAux(1, 0, 0, 8, .075, .075, 'rand', 10000, 1000, 2, 1000, 30)
def ComputeInOutPOCorrParallelAux(kappa, p, gamma, nPhis, mExt, mExtOne, rewireType, N, K, nPop, T, trNo):
    out = np.nan
    try:
	tcIn = GetInputTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa, 'net')
	tcOut = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
	theta = np.linspace(0, 180, nPhis, endpoint = False)
	poIn = POofPopulation(tcIn, theta, IF_IN_RANGE = True) * np.pi / 180
	poOut = POofPopulation(tcOut, theta, IF_IN_RANGE = True) * np.pi / 180
	print 'computing ccc... ', trNo
	sys.stdout.flush()
	# plt.plot(poIn[:N], poOut[:N], '.')
	# plt.show()
	# ipdb.set_trace()
	out = CircularCorrCoeff(poIn[:N], poOut[:N])
	print 'done', trNo
    except IOError:
	print ''
	return out
    return out
    
def PlotInOutPOCorrParallel(nRewireSteps, kappa = 1, p = 0, gamma = 0, nPhis = 8, mExt = .075, mExtOne = .075, rewireType = 'rand', trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_COMPUTE = True):
    outE = []; outI = []
    validTr = []
    if rewireType == 'rand':
	filename = './data/twopop/InOutCorr_p%sg%sk%s'%(int(100 * p), int(10 * gamma), int(10 * kappa))
    else:
	filename = './data/twopop/InOutCorr_p%sg%sk%s'%(int(100 * p), int(10 * gamma), int(10 * kappa)) + '_'  + rewireType 
	
    if IF_COMPUTE:
        pool = Pool(nRewireSteps)	
        func = partial(ComputeInOutPOCorrParallelAux, kappa, p, gamma, nPhis, mExt, mExtOne, rewireType, N, K, nPop, T)
	outE = pool.map(func, range(nRewireSteps))
	# ipdb.set_trace()
        np.save(filename, [validTr, outE])
	pool.close()
    else:
	filename += '.npy'
	outE = np.asarray(np.load(filename)[1], dtype = float)
	# ipdb.set_trace()	
	validTr =  np.arange(nRewireSteps)
	validTr = validTr[~np.isnan(outE)]
	plt.figure()
	plt.plot(validTr, outE, 'ko-')
	plt.xlabel('rewire step')
	# plt.ylabel(r'$\langle \cos 2 [ \phi_{j}^{po} - \theta_{j, IN}^{po} ]   \rangle_j$')
	plt.ylabel('CCC')
	plt.plot(outI, 'ro-')
	filepath = './figs/twopop/InOutCorr_p%sg%s'%(int(100 * p), int(10 * gamma)) + '_' + rewireType
	ProcessFigure(plt.gcf(), filepath, True)

def ProcessFigure(figHdl, filepath, IF_SAVE, IF_XTICK_INT = False, figFormat = 'eps', paperSize = [4, 3], titleSize = 10, axPosition = [0.25, 0.25, .65, .65], tickFontsize = 10, labelFontsize = 12, nDecimalsX = 1, nDecimalsY = 1):
    # FixAxisLimits(figHdl)
    # FixAxisLimits(plt.gcf(), IF_XTICK_INT, nDecimalsX, nDecimalsY)
    Print2Pdf(plt.gcf(), filepath, paperSize, figFormat=figFormat, labelFontsize = labelFontsize, tickFontsize=tickFontsize, titleSize = titleSize, IF_ADJUST_POSITION = True, axPosition = axPosition)
    plt.show()

def FixAxisLimits(fig, IF_XTICK_INT = False, nDecimalsX = 1, nDecimalsY = 1):
    ax = fig.axes[0]
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.' + '%s'%(int(nDecimalsX)) + 'f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.' + '%s'%(int(nDecimalsY)) + 'f'))
    xmiddle = 0.5 * (xmin + xmax)
    xticks = [xmin, xmiddle, xmax]
    if IF_XTICK_INT:
	if xmiddle != int(xmiddle):
	    xticks = [xmin, xmax]
	ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.set_xticks(xticks)
    ax.set_yticks([ymin, 0.5 *(ymin + ymax), ymax])
    plt.draw()

def press(event):
    print('press', event.key)
    sys.stdout.flush()
    if event.key == 'w':
	neuronid = int(plt.gca().get_title())
	filename = "./figs/twopop/rewire/neuron%s"%(neuronid)
	paperSize = [4, 3]
        Print2Pdf(plt.gcf(),  filename,  paperSize, figFormat='png', labelFontsize = 10, tickFontsize=8, titleSize = 10.0, IF_ADJUST_POSITION = True, axPosition = [0.14, 0.14, .7, .7])
	fig.canvas.draw()

def FollowNeurons(neuronsList, p, gamma, nPhis, mExt, mExtOne, rewireType, trList = [0], N = 10000, K = 1000, nPop = 2, T = 1000, IF_LEGEND = True):
        tc = [[] for x in xrange(len(trList))]
	preferredOri = [[] for x in xrange(len(trList))]
	theta = np.linspace(0, 180, nPhis, endpoint = False)
	for kk, trNo in enumerate(trList):
	    tmp = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T)
	    tc[kk].append(tmp)
	    preferredOri[kk].append(POofPopulation(tc[kk][0][:N], IF_IN_RANGE = True))
        # ipdb.set_trace()
        fig, ax = plt.subplots()
	plt.ion()
	fig.canvas.mpl_connect('key_press_event', press)
	for neuronIdx in neuronsList:
	    for kk, trNo in enumerate(trList):
		osi = OSI(tc[kk][0][neuronIdx, :], theta)
		plt.plot(theta, tc[kk][0][neuronIdx, :], 's-', label = 'STP:%s OSI:%.4s PO:%.5s'%(trNo, osi, preferredOri[kk][0][neuronIdx]))
		plt.gca().set_position([0.15, 0.15, .65, .65])
		paperSize = [4, 3]
                FixAxisLimits(plt.gcf())
		plt.title('%s'%(neuronIdx))
		plt.xlabel('Stimulus (deg)')
		plt.ylabel(r'$m(\phi)$')
		plt.show()
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
	plt.title('%s'%(neuronIdx))
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

# def GetPOofPop(p, gamma, mExt, mExtOne, rewireType = 'rand', nPhis = 8, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_IN_RANGE = True):
#     nNeurons = N
#     thetas= np.linspace(0, np.pi, nNeurons, endpoint = False)
#     # ipdb.set_trace()
#     tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T)
#     prefferedOri = POofPopulation(tc[:N], IF_IN_RANGE = True) * np.pi / 180.0
#     return prefferedOri

def GetPOofPop(p, gamma, mExt, mExtOne, rewireType, nPhis = 8, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_IN_RANGE = True, kappa = 1):
    nNeurons = N
    thetas= np.linspace(0, np.pi, nNeurons, endpoint = False)
    tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
    prefferedOri = POofPopulation(tc[:N], IF_IN_RANGE = True) * np.pi / 180.0
    return prefferedOri

def GetPOofPopAllNeurons(p, gamma, mExt, mExtOne, rewireType, nPhis = 8, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_IN_RANGE = True, kappa = 1):
    nNeurons = N
    thetas= np.linspace(0, np.pi, nNeurons, endpoint = False)
    tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
    print tc.shape
    prefferedOri = POofPopulation(tc, IF_IN_RANGE = True) * np.pi / 180.0

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
    atThetas = np.linspace(0, 180, nThetas, endpoint = False)
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

def TrackPOofPop(trList, p, gamma, mExt, mExtOne, rewireType = 'rand', nPhis = 8, N = 10000, K = 1000, nPop = 2, T = 1000, IF_IN_RANGE = True, kappa = 1):
    po = []
    z = []
    plt.figure()
    for trNo in trList:
	po.append(GetPOofPop(p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N, K, nPop, T, IF_IN_RANGE, kappa))
    for i in range(1, len(trList)):
	z.append(np.nanmean(np.cos(2 * (po[i - 1][:N] - po[i][:N]))))
    plt.plot(z, 'k.-')
    FixAxisLimits(plt.gcf())
    plt.gca().set_position([0.25, 0.25, .65, .65])
    plt.xlim(0, trList[-1])
    paperSize = [4, 3]
    FixAxisLimits(plt.gcf())
    plt.title('PO stability')
    plt.xlabel('rewiring step (n)')
    plt.ylabel(r'$\langle \cos 2 [ \phi_{i, n}^{po} - \phi_{i, n + 1}^{po} ]   \rangle_i$')
    filename = ''
    filename = filename + "p%sg%s"%(p, gamma)
    Print2Pdf(plt.gcf(),  "./figs/twopop/PO_vs_n"+filename,  paperSize, figFormat='png', labelFontsize = 12, tickFontsize=8, titleSize = 10.0)
    plt.show()
    return z    

def PltOSIHist(p, gamma, nPhis=8, mExt=.075, mExtOne=.075, rewireType = 'rand', trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_NEW_FIG = True, color = 'k', kappa = 1, legendTxt = ''):
    NE = N
    NI = N
    tc = np.zeros((NE + NI, nPhis))
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    if IF_NEW_FIG:
        fg0, ax0 = plt.subplots()
        fg1, ax1 = plt.subplots()
    else:
        fg0 = plt.figure(1)
        fg1 = plt.figure(2)
    for i, iPhi in enumerate(phis):
	print i, iPhi,
	try:
	    if i == 0:
		print 'loading from fldr: ',	    
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = True, kappa = kappa)
                # ipdb.set_trace()
	    else:
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = False, kappa = kappa)
	    if(len(fr) == 1):
		if(np.isnan(fr)):
		    print 'file not found!'
	    tc[:, i] = fr
            print ''
	except IOError:
	    print 'file not found!'
    osi = OSIOfPop(tc[:NE, :], phis)
    osiI = OSIOfPop(tc[NE:, :], phis)
    print "K = ", K, ", osi simulation: ", np.nanmean(osi)
    plt.figure(fg0.number)
    # plt.xlabel(r"$\mathrm{OSI} \,\,\,\,  (m_{E, i}^{(1)})$")
    plt.xlabel('OSI', fontsize = 12)    
    plt.ylabel('Density', fontsize = 12)
    # plt.hist(osi[~np.isnan(osi)], 27, normed = 1, histtype = 'step', label = r'$STP: %s$'%(trNo, ) + legendTxt, color = color, lw = 1)
    plt.hist(osi[~np.isnan(osi)], 27, normed = 1, histtype = 'step', label = legendTxt, color = color, lw = 1)

    # plt.hist(osiI[~np.isnan(osiI)], 27, normed = 1, histtype = 'step', label = legendTxt, color = 'r', lw = 1)

    
    # plt.hist(osi[~np.isnan(osi)], 27, normed = 1, histtype = 'step', label = r'$p = %s,\,\gamma = %s$'%(p, gamma, ), color = color, lw = 1)
    plt.xlim(0, 1)
    plt.gca().set_xticks([0, 0.5, 1])
    _, ymax = plt.ylim()
    plt.gca().set_yticks([0, np.ceil(ymax)])    
    # plt.title(r'$N = %s,\, K = %s,\, m_0^{(0)} = %s,\, m_0^{(1)} = %s$'%(NE, K, mExt, mExtOne))
    # plt.title(r'$m_0^{(0)} = %s, \,m_0^{(1)} = %s $'%(mExt, mExtOne))
    # plt.title(r'$E, p = %s, \, \gamma = %s $'%(p, gamma))
    _, ymax = plt.ylim()
    IF_PLOT_VLINE = False
    if IF_PLOT_VLINE:
	plt.vlines(np.nanmean(osi), 0, ymax, lw = 1, color = color)
    ax = plt.gca()
    # ann = ax.annotate('', xy=(np.nanmean(osi), 0), xytext=(np.nanmean(osi), 1.5),arrowprops=dict(facecolor='black', arrowstyle = 'simple', color = color))
    print "mean OSI, E = ", np.nanmean(osi)
    print "mean OSI, I= ", np.nanmean(osiI)
    # plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 8})
    filename = './figs/twopop/' + 'OSIhist_E_vs_kappa_p%s_gamma%s_m0%s'%(p, gamma, int(mExt*1e3))
    paperSize = [2.5, 2]
    axPosition=[.26, .24, .65, .65]
    ProcessFigure(plt.gcf(), filename, 1, paperSize = paperSize, axPosition = axPosition, titleSize=10, nDecimalsX=1, nDecimalsY=1, figFormat='svg', labelFontsize = 8, tickFontsize = 6)
    plt.figure(fg1.number)    
    # plt.title(r'$I, p = %s, \, \gamma = %s $'%(p, gamma))
    plt.hist(osiI[~np.isnan(osiI)], 27, normed = 1, histtype = 'step', label = legendTxt, color = color, lw = 1)
    ax = plt.gca()
    # ann = ax.annotate('', xy=(np.nanmean(osiI), 0), xytext=(np.nanmean(osiI), 1.5),arrowprops=dict(facecolor='black', arrowstyle = 'simple', color = color))
    plt.xlabel('OSI', fontsize = 12)    
    plt.ylabel('Density', fontsize = 12)
    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 6})
    filename = './figs/twopop/' + 'OSIhist_I_vs_kappa_p%s_gamma%s_m0%s'%(p, gamma, int(mExt*1e3))
    paperSize = [2.5, 2]
    axPosition=[.26, .24, .65, .65]
    plt.xlim(0, 1)
    plt.gca().set_xticks([0, 0.5, 1])
    _, ymax = plt.ylim()
    plt.gca().set_yticks([0, np.ceil(ymax)])    
    ProcessFigure(plt.gcf(), filename, 1, paperSize = paperSize, axPosition = axPosition, titleSize=10, nDecimalsX=1, nDecimalsY=1, figFormat='svg', labelFontsize = 8, tickFontsize = 6)
    return osi, osiI

def CompareMeanOSIHist(pList, gList, nPhis = 8, mExt = 0.075, mExtOneList  = [0.075], trList = [0], rewireType = 'cntrl', N = 10000, KList = [1000], nPop = 2, T = 1000, IF_NEW_FIG = True, clrCntr = 0, filename = '', IF_LEGEND = True, legendTxt = '', kappa = 1, color = '', neuronType = 'E'):
    # colors = [plt.cm.Dark2(i) for i in np.linspace(0, 1, 1 + clrCntr + len(pList) * len(gList) * len(mExtOneList) * len(trList), endpoint = False)]
    colors = [plt.cm.Dark2(i) for i in np.linspace(0, 1, clrCntr + len(pList) * len(gList) * len(mExtOneList) * len(KList), endpoint = False)]
    meanOSI = []
    meanOSI_I = []
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    validTrList = []
    legendTxtTag = legendTxt
    if IF_NEW_FIG:
	plt.figure()
    for mExtOne in mExtOneList:
	for kk, K in enumerate(KList):
	    legendTxt = ', K=%s'%(K) + legendTxtTag
	    validTrList = []
	    meanOSI = []
	    meanOSI_I = []
	    counter = 0
	    for p in pList:
		for gamma in gList:
		    for trNo in trList:
			try:
			    print trNo
                            tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
			    osi = OSIOfPop(tc[:N, :], phis)
			    osiI = OSIOfPop(tc[N:, :], phis)
			    meanOSI.append(np.nanmean(osi))
			    meanOSI_I.append(np.nanmean(osiI))
			    if ~np.isnan(meanOSI[counter]):
				validTrList.append(trNo)
			    clrCntr += 1
			    counter += 1
			except IOError:
			    print "p = ", p, " gamma = ", gamma, " trial# ", trNo, " file not found"
	    # ipdb.set_trace()
	    validTrList = np.asarray(validTrList, dtype = int)
	    meanOSI = np.array(meanOSI)
	    meanOSI_I = np.array(meanOSI_I)
	    pcolor = color
	    if color == '':
	    	pcolor = colors[kk]
	    if neuronType == 'E':
		plt.plot(validTrList, meanOSI[validTrList], 'o-', color = pcolor, label = r'$\kappa=%s$'%(kappa) + legendTxt, markeredgecolor = pcolor)
	    else:
		plt.plot(validTrList, meanOSI_I[validTrList], 'o-', color = pcolor, label = r'$\kappa=%s$'%(kappa) + legendTxt, markeredgecolor = pcolor)

	    
 
    print '--'*26
    osiLast = meanOSI[-1]
    osilastCnt = -1
    osiLastI = meanOSI_I[-1]
    while np.isnan(osiLast):
	osilastCnt -= 1
	print osilastCnt
	osiLast = meanOSI[osilastCnt]	
	osiLastI = meanOSI_I[osilastCnt]
    print 'pc change in mean OSI = ', 100 * (osiLast - meanOSI[0]) / meanOSI[0]
    print '--'*26
    plt.gca().set_position([0.15, 0.15, .65, .65])
    # ipdb.set_trace()
    # if neuronType == 'E':
    # 	plt.plot(validTrList, meanOSI, '*-', color = color, label = r'$\kappa=%s$'%(kappa) + legendTxt)
    # 	plt.title('E')
    # else:
    # 	plt.plot(validTrList, meanOSI_I, '*-', color = color, label = r'$\kappa=%s$'%(kappa) + legendTxt)
    # 	plt.title('I')
    plt.ylim(0, .5)
    xmax = len(trList) - 1
    # if np.mod(xmax, 2) != 0:
    # 	xmax += 1
    plt.xlim(0, xmax)
    plt.xlabel('rewiring step')
    plt.ylabel(r'$\langle OSI \rangle$')
    plt.gca().set_position([0.25, 0.25, .65, .65])
    filename = filename + "p%sg%sk%s_K%s_"%(p, gamma, kappa, K) + neuronType
    filepath = "./figs/twopop/compareOSI_mean_"+ rewireType + '_' +  filename 
    paperSize = [4, 3]
    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 8})
    plt.grid('on')
#    ipdb.set_trace()
    print 'saving as: ', filepath
    ProcessFigure(plt.gcf(), filepath, IF_SAVE = 1, IF_XTICK_INT = True, figFormat = 'eps')
    
    #Print2Pdf(plt.gcf(),  "./figs/twopop/compareOSI_mean_"+filename,  paperSize, figFormat='pdf', labelFontsize = 10, tickFontsize=8, titleSize = 10.0)
    plt.show()
    print meanOSI


def CompareOSIHistFig_sup1(pList, gList, nPhis = 8, mExt= .075, mExtOneList = [.075], trList= [0], rewireType = 'rand', N = 10000, KList = [1000], nPop = 2, T = 1000, IF_NEW_FIG = True, clrCntr = 0, filename = '', IF_LEGEND = True, legendTxt = '', kappa = 1, color = ''):
    if IF_NEW_FIG:
	plt.figure()
    if color == '':
	colors = [plt.cm.Dark2(i) for i in np.linspace(0, 1, 1 + clrCntr + len(pList) * len(gList) * len(mExtOneList) * len(trList), endpoint = False)]
    else:
	pcolor = color
    meanOSI = []
    meanOSI_I = []


    K=1000
    kappa = 0
    trNo = 100
    
    tmposi, tmposiI = PltOSIHist(p=0, gamma=0, nPhis=8, mExt=0.075, mExtOne=0.075, rewireType='rand', trNo = trNo, IF_NEW_FIG = False, color = 'k', T=T, K=K, kappa = kappa, legendTxt = legendTxt, N=N)

    kappa = 8
    trNo = 2000
    tmposi, tmposiI = PltOSIHist(p=0, gamma=0, nPhis=8, mExt=0.075, mExtOne=0.075, rewireType='rand', trNo = trNo, IF_NEW_FIG = False, color = 'g', T=T, K=K, kappa = kappa, legendTxt = legendTxt, N=N)
    

    ipdb.set_trace()



    
    # for mExtOne in mExtOneList:
    #     for trNo in trList:
    #         for p in pList:
    #     	for gamma in gList:
    #     	    for K in KList:
    #     		try:
    #     		    print trNo
    #     		    if color == '':
    #     			tmposi, tmposiI = PltOSIHist(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo = trNo, IF_NEW_FIG = False, color = colors[clrCntr], T=T, K=K, kappa = kappa, legendTxt = legendTxt, N=N)
    #     		    else:
    #     			tmposi, tmposiI = PltOSIHist(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo = trNo, IF_NEW_FIG = False, color = pcolor, T=T, K=K, kappa = kappa, legendTxt = legendTxt, N=N)
    #     		    meanOSI.append(np.nanmean(tmposi))
    #     		    meanOSI_I.append(np.nanmean(tmposiI))
    #     		    clrCntr += 1
    #     		except IOError:
    #     		    print "p = ", p, " gamma = ", gamma, " trial# ", trNo, " file not found"
    # plt.gca().legend(bbox_to_anchor = (1.1, 1.5))
    if IF_LEGEND:
	plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 8})
    # print '--'*26
    # osiLast = meanOSI[-1]
    # osilastCnt = -1
    # osiLastI = meanOSI[-1]
    # while np.isnan(osiLast):
    #     osilastCnt -= 1
    #     print osilastCnt
    #     osiLast = meanOSI[osilastCnt]	
    #     osiLastI = meanOSI_I[osilastCnt]
    # print 'pc change in mean OSI = ', 100 * (osiLast - meanOSI[0]) / meanOSI[0]
    # # print 'pc change in mean OSI, I = ', 100 * (osiLastI - meanOSI_I[0]) / meanOSI_I[0]
    # print '--'*26
    plt.gca().set_position([0.15, 0.15, .65, .65])
    if nPop == 2:
	# plt.savefig("./figs/twopop/compareOSI_.png")
	# plt.savefig("./figs/twopop/compareOSI_"+filename + '.png')
	paperSize = [4, 3]
	# ipdb.set_trace()
        p=0
        gamma=0
	filename = filename + "p%sg%sk%s"%(p, gamma, kappa)
        plt.figure(1)
        plt.gca().set_xticklabels(['0', '0.5', '1'])
        plt.draw()        
        Print2Pdf(plt.gcf(),  "./figs/twopop/compareOSI_" + rewireType + '_' + filename,  paperSize, figFormat='svg', labelFontsize = 10, tickFontsize=8, titleSize = 10.0, IF_ADJUST_POSITION = True, axPosition = [0.14, 0.14, .7, .7])
    # plt.figure()
    # plt.plot(trList, meanOSI, 'k*-')
    plt.xlabel('rewiring step')
    plt.ylabel(r'$\langle OSI \rangle$')
    plt.gca().set_position([0.25, 0.25, .65, .65])
    # plt.xlim([0, ])
    filename = filename + "p%sg%s"%(p, gamma)
    plt.figure(2)
    plt.gca().set_xticklabels(['0', '0.5', '1'])
    plt.draw()
    Print2Pdf(plt.gcf(),  "./figs/twopop/compareOSI_mean_"+ rewireType + '_' + filename,  paperSize, figFormat='svg', labelFontsize = 10, tickFontsize=8, titleSize = 10.0)
    plt.show()
    print meanOSI

    

def CompareOSIHist(pList, gList, nPhis = 8, mExt= .075, mExtOneList = [.075], trList= [0], rewireType = 'rand', N = 10000, KList = [1000], nPop = 2, T = 1000, IF_NEW_FIG = True, clrCntr = 0, filename = '', IF_LEGEND = True, legendTxt = '', kappa = 1, color = ''):
    if IF_NEW_FIG:
	plt.figure()
    if color == '':
	colors = [plt.cm.Dark2(i) for i in np.linspace(0, 1, 1 + clrCntr + len(pList) * len(gList) * len(mExtOneList) * len(trList), endpoint = False)]
    else:
	pcolor = color
    meanOSI = []
    meanOSI_I = []
    for mExtOne in mExtOneList:
	for trNo in trList:
	    for p in pList:
		for gamma in gList:
		    for K in KList:
			try:
			    print trNo
			    if color == '':
				tmposi, tmposiI = PltOSIHist(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo = trNo, IF_NEW_FIG = False, color = colors[clrCntr], T=T, K=K, kappa = kappa, legendTxt = legendTxt, N=N)
			    else:
				tmposi, tmposiI = PltOSIHist(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo = trNo, IF_NEW_FIG = False, color = pcolor, T=T, K=K, kappa = kappa, legendTxt = legendTxt, N=N)
			    meanOSI.append(np.nanmean(tmposi))
			    meanOSI_I.append(np.nanmean(tmposiI))
			    clrCntr += 1
			except IOError:
			    print "p = ", p, " gamma = ", gamma, " trial# ", trNo, " file not found"
    # plt.gca().legend(bbox_to_anchor = (1.1, 1.5))
    if IF_LEGEND:
	plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 8})
    print '--'*26
    osiLast = meanOSI[-1]
    osilastCnt = -1
    osiLastI = meanOSI[-1]
    while np.isnan(osiLast):
	osilastCnt -= 1
	print osilastCnt
	osiLast = meanOSI[osilastCnt]	
	osiLastI = meanOSI_I[osilastCnt]
    print 'pc change in mean OSI = ', 100 * (osiLast - meanOSI[0]) / meanOSI[0]
    # print 'pc change in mean OSI, I = ', 100 * (osiLastI - meanOSI_I[0]) / meanOSI_I[0]
    print '--'*26
    plt.gca().set_position([0.15, 0.15, .65, .65])
    if nPop == 2:
	# plt.savefig("./figs/twopop/compareOSI_.png")
	# plt.savefig("./figs/twopop/compareOSI_"+filename + '.png')
	paperSize = [4, 3]
	# ipdb.set_trace()
	filename = filename + "p%sg%sk%s"%(p, gamma, kappa)
        Print2Pdf(plt.gcf(),  "./figs/twopop/compareOSI_" + rewireType + '_' + filename,  paperSize, figFormat='eps', labelFontsize = 10, tickFontsize=8, titleSize = 10.0, IF_ADJUST_POSITION = True, axPosition = [0.14, 0.14, .7, .7])
    # plt.figure()
    # plt.plot(trList, meanOSI, 'k*-')
    plt.xlabel('rewiring step')
    plt.ylabel(r'$\langle OSI \rangle$')
    plt.gca().set_position([0.25, 0.25, .65, .65])
    # plt.xlim([0, ])
    filename = filename + "p%sg%s"%(p, gamma)
    Print2Pdf(plt.gcf(),  "./figs/twopop/compareOSI_mean_"+ rewireType + '_' + filename,  paperSize, figFormat='png', labelFontsize = 10, tickFontsize=8, titleSize = 10.0)
    plt.show()
    print meanOSI

def MeanRateVsRewireStep(kappa = 1, p = 0, gamma = 0, nPhis = 8, mExt = 0.075, mExtOne = .075, trList = range(30), rewireType = 'rand', N = 10000, K = 1000, nPop = 2, T = 1000, IF_NEW_FIG = True, clrCntr = 0, filename = '', IF_LEGEND = True, legendTxt = '', paperSize = [1.6, 1.2], axPosition = [0.3, 0.3, .5, .5]):
    if IF_NEW_FIG:
	plt.figure()
#    colors = [plt.cm.Dark2(i) for i in np.linspace(0, 1, 1 + clrCntr + len(pList) * len(gList) * len(mExtOneList) * len(trList), endpoint = False)]
    meanRateE = []; meanRateI = []
    for trNo in trList:
	print trNo
	try:
	    print trNo
	    tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
	    meanRateE.append(tc[:N].mean())
	    meanRateI.append(tc[N:].mean())
	except IOError:
	    print ''
    plt.plot(trList, meanRateE, 'k.-', label = 'E')
    plt.plot(trList, meanRateI, 'r.-', label = 'I')
    plt.legend(loc = 0, frameon = False, numpoints = 1, ncol = 2, prop = {'size': 8})        
    plt.xlabel('rewire step')
    # plt.ylabel(r'$\left[ \langle m_i(\phi) \rangle_i \right]_{\phi}$')
    plt.ylabel('mean rate')
    plt.title(r'$\kappa = %s$'%(kappa))
    plt.grid('on')
    plt.ylim(0, .5)
    plt.xlim(0, 6)
    filename = './figs/twopop/meanrate_vs_rewire_p%sg%sk%s'%(int(100 * p), int(10 * gamma), int(10 * kappa))
    ProcessFigure(plt.gcf(), filename, True, True, figFormat='eps', paperSize = paperSize, axPosition = axPosition, tickFontsize = 6, labelFontsize = 8)
    print filename
    ipdb.set_trace()

    # for mExtOne in mExtOneList:
    # 	    for p in pList:
    # 		for gamma in gList:
    # 		    for K in KList:

def ComputeFFInOutPOCorrParallelAux(kappa, p, gamma, nPhis, mExt, mExtOne, rewireType, N, K, nPop, T, NFF, JE0, JI0, IF_PLOT, trNo):
    out = np.nan
    try:
	theta = np.linspace(0, 180, nPhis, endpoint = False)
	uFF = rw.ComputeFFInput(nPhis, p, gamma, kappa, mExt, mExtOne, trNo, N, K, nPop, NFF, JE0, JI0, 0.2, rewireType, T) #nPhis, p, gamma, kappa, mExt, mExtOne, trNo, nPop, NFF, JE0, JI0, cFF, rewireType, T)
	# tcIn = GetInputTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa, 'FF')
	# ipdb.set_trace()
        # ipdb.set_trace()	
        poFF = POofPopulation(uFF, theta, IF_IN_RANGE = True) * np.pi / 180
	tcOut = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
	poOut = POofPopulation(tcOut, theta, IF_IN_RANGE = True) * np.pi / 180
	if IF_PLOT:
	    plt.plot(poFF[:N], poOut[:N], '.')
	    plt.show()
	    ipdb.set_trace()
	print 'computing ccc... ', trNo
	sys.stdout.flush()
	out = CircularCorrCoeff(poFF[:N], poOut[:N])
	print 'done', trNo, 'CCC=', out
    except IOError:
	print ''
	return out
    return out

def PlotFFInOutPOCorrParallel(nRewireSteps, kappa = 1, p = 0, gamma = 0, nPhis = 8, mExt = .075, mExtOne = .075, rewireType = 'rand', trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_COMPUTE = True, JE0 = 2.0, JI0 = 1.0, NFF = 10000):
    outE = []; outI = []
    validTr = []
    if rewireType == 'rand':
	filename = './data/twopop/FFInOutPOCorr_p%sg%sk%s'%(int(100 * p), int(10 * gamma), int(10 * kappa))
    else:
	filename = './data/twopop/FFInOutPOCorr_p%sg%sk%s'%(int(100 * p), int(10 * gamma), int(10 * kappa)) + '_' + rewireType
    if IF_COMPUTE:
        pool = Pool(nRewireSteps)	
        func = partial(ComputeFFInOutPOCorrParallelAux, kappa, p, gamma, nPhis, mExt, mExtOne, rewireType, N, K, nPop, T, NFF, JE0, JI0, False)
	outE = pool.map(func, range(nRewireSteps))
        np.save(filename, [validTr, outE])
	pool.close()
    else:
	filename += '.npy'
	outE = np.asarray(np.load(filename)[1], dtype = float)
	validTr =  np.arange(nRewireSteps)
	validTr = validTr[~np.isnan(outE)]
	plt.figure()
	plt.plot(validTr, outE, 'ko-')
	plt.xlabel('rewire step')
	plt.title('FF PO vs Out PO')
	# plt.ylabel(r'$\langle \cos 2 [ \phi_{j}^{po} - \theta_{j, IN}^{po} ]   \rangle_j$')
	plt.ylabel('CCC')
	plt.plot(outI, 'ro-')
	filepath = './figs/twopop/FFInOutPOCorr_p%sg%s'%(int(100 * p), int(10 * gamma)) + '_' + rewireType
	ProcessFigure(plt.gcf(), filepath, True)

def ComputePOCorrs(kappa, nRewireSteps, rewireType):
    cccNet = PlotInOutPOCorrParallel(nRewireSteps, kappa=kappa, rewireType = rewireType)
    cccFF = PlotFFInOutPOCorrParallel(nRewireSteps, IF_COMPUTE = 1, kappa=kappa, rewireType = rewireType)
    cccRec = PlotRecInOutPOCorrParallel(nRewireSteps, IF_COMPUTE = 1, kappa=kappa, rewireType = rewireType)
    
def PlotInOutPOCorrCmpnts(nRewireSteps, color = 'k', kappa = 1, p = 0, gamma = 0, nPhis = 8, mExt = .075, mExtOne = .075, rewireType = 'rand', trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, JE0 = 2.0, JI0 = 1.0, NFF = 10000, IF_NEW_FIG = True, startIdx = 0):
    outE = []; outI = []
    validTr = []
    if rewireType == 'rand':
        filename = './data/twopop/FFInOutPOCorr_p%sg%sk%s'%(int(100 * p), int(10 * gamma), int(10 * kappa))
    if rewireType == 'twosteps':
        filename = './data/twopop/FFInOutPOCorr_p%sg%sk%s_twosteps'%(int(100 * p), int(10 * gamma), int(10 * kappa))	
    
    filename += '.npy'
    outE = np.asarray(np.load(filename)[1], dtype = float)
    outE = outE[:nRewireSteps]
    validTr =  np.arange(nRewireSteps, dtype = int)
    validTr = validTr[~np.isnan(outE)]
    if IF_NEW_FIG:
	plt.figure()
    plt.plot(validTr, outE, 'o-', label = r'$\mathrm{FF}, \kappa=%s$'%(kappa), markersize = 5, color=color, markeredgecolor=color)
    plt.xlabel('rewire step')
    plt.ylabel('Circ Corr Coeff')
    plt.plot(outI, 'ro-')
    if rewireType == 'rand':
	filename = './data/twopop/RecInOutPOCorr_p%sg%sk%s'%(int(100 * p), int(10 * gamma), int(10 * kappa))
    if rewireType == 'twosteps':
	filename = './data/twopop/RecInOutPOCorr_p%sg%sk%s_twosteps'%(int(100 * p), int(10 * gamma), int(10 * kappa))
	
    filename += '.npy'
    recCurPOCorr = np.asarray(np.load(filename)[1], dtype = float)[:nRewireSteps]
    plt.plot(recCurPOCorr, 'k*-', label = r'$\mathrm{rec}, \kappa=%s$'%(kappa), markersize = 5, color=color, markeredgecolor=color)
    if rewireType == 'rand':
	filename = './data/twopop/InOutCorr_p%sg%sk%s'%(int(100 * p), int(10 * gamma), int(10 * kappa))
    if rewireType == 'twosteps':
	filename = './data/twopop/InOutCorr_p%sg%sk%s_twosteps'%(int(100 * p), int(10 * gamma), int(10 * kappa))
	
    filename += '.npy'
    netCurPOCorr = np.asarray(np.load(filename)[1], dtype = float)[:nRewireSteps]
    validTr =  np.arange(nRewireSteps)
    validTr = validTr[~np.isnan(netCurPOCorr[:nRewireSteps])]
    plt.plot(netCurPOCorr, 'ks-', label =  r'$\mathrm{net}, \kappa=%s$'%(kappa), markersize = 5, color=color, markeredgecolor=color)    
    plt.legend(loc = 0, frameon = False, numpoints = 1, ncol = 2, prop = {'size': 8})    
    filepath = './figs/twopop/FF_nd_InOutPOCorr_p%sg%sk%s_'%(int(100 * p), int(10 * gamma), int(10*kappa)) + rewireType
    print 'saving figure as:', filepath
    # plt.gca().set_xticks([0, 15, 30])
    plt.ylim(0, 1)
    plt.gca().set_yticks([0, 0.5, 1])
    plt.title('CCC: PO out vs PO in components')
    ProcessFigure(plt.gcf(), filepath, True, True)

def PlotRecInOutPOCorrParallel(nRewireSteps, kappa = 1, p = 0, gamma = 0, nPhis = 8, mExt = .075, mExtOne = .075, rewireType = 'rand', trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_COMPUTE = True, JE0 = 2.0, JI0 = 1.0, NFF = 10000):
    outE = []; outI = []
    validTr = []
    if rewireType == 'rand':
	filename = './data/twopop/RecInOutPOCorr_p%sg%sk%s'%(int(100 * p), int(10 * gamma), int(10 * kappa))
    else:
	filename = './data/twopop/RecInOutPOCorr_p%sg%sk%s'%(int(100 * p), int(10 * gamma), int(10 * kappa)) + '_' + rewireType

    if IF_COMPUTE:
        pool = Pool(nRewireSteps)	
        func = partial(ComputeRecInOutPOCorrParallelAux, kappa, p, gamma, nPhis, mExt, mExtOne, rewireType, N, K, nPop, T, NFF, JE0, JI0, False)
	outE = pool.map(func, range(nRewireSteps))
        np.save(filename, [validTr, outE])
	pool.close()
    else:
	filename += '.npy'
	outE = np.asarray(np.load(filename)[1], dtype = float)
	validTr =  np.arange(nRewireSteps)
	validTr = validTr[~np.isnan(outE)]
	plt.figure()
	plt.plot(validTr, outE, 'ko-')
	plt.xlabel('rewire step')
	plt.ylabel('CCC')
	plt.plot(outI, 'ro-')
	filepath = './figs/twopop/Rec_FF_nd_Net_InOutPOCorr_p%sg%s'%(int(100 * p), int(10 * gamma)) + '_' + rewireType
	ProcessFigure(plt.gcf(), filepath, True)

def ComputeRecInOutPOCorrParallelAux(kappa, p, gamma, nPhis, mExt, mExtOne, rewireType, N, K, nPop, T, NFF, JE0, JI0, IF_PLOT, trNo):
    out = np.nan
    try:
	theta = np.linspace(0, 180, nPhis, endpoint = False)
	uFF = rw.ComputeFFInput(nPhis, p, gamma, kappa, mExt, mExtOne, trNo, N, K, nPop, NFF, JE0, JI0, 0.2, rewireType, T) #nPhis, p, gamma, kappa, mExt, mExtOne, trNo, nPop, NFF, JE0, JI0, cFF, rewireType, T)
	tcIn = GetInputTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
	# ipdb.set_trace()
	uRec = tcIn - uFF
        # ipdb.set_trace()	
        poRec = POofPopulation(uRec, theta, IF_IN_RANGE = True) * np.pi / 180
	tcOut = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
	poOut = POofPopulation(tcOut, theta, IF_IN_RANGE = True) * np.pi / 180
	if IF_PLOT:
	    plt.plot(poFF[:N], poOut[:N], '.')
	    plt.show()
	    ipdb.set_trace()
	print 'computing ccc... ', trNo
	sys.stdout.flush()
	out = CircularCorrCoeff(poRec[:N], poOut[:N])
	print 'done', trNo, 'CCC=', out
    except IOError:
	print ''
	return out
    return out

def LocalPOCorrFunc(po, idx, nNeighbours, N):
    neighbourIdx = np.mod(np.arange(idx - nNeighbours, idx + nNeighbours + 1), N)
    # ipdb.set_trace()
    # return np.nanmean(CircularCorrCoeff(po[idx],  po[neighbourIdx[neighbourIdx != idx]]))
    return np.nanmean(np.cos(2 * (po[idx] - po[neighbourIdx[neighbourIdx != idx]])))
    
def LocalPOCorr(kappa=1, p=0, gamma=0, nPhis=8, mExt=0.075, mExtOne=0.075, rewireType='rand', N=10000, K=1000, nPop=2, T=1000, trNo=0, nNeighbours = 100):
    theta = np.linspace(0, 180, nPhis, endpoint = False)
    tcOut = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
    poOut = POofPopulation(tcOut, theta, IF_IN_RANGE = True) * np.pi / 180
    poOut = poOut[:N]
    localCorr = []
    [localCorr.append(LocalPOCorrFunc(poOut, idx, nNeighbours, N)) for idx in xrange(N)]
    return np.asarray(localCorr)

def LoadM1vsT(p = 0, gamma = 0, phi = 0, trNo = 0, mExt = 0.075, mExtOne = 0.075, K = 1000, NE = 10000, T = 1000, nPop = 2, rewireType = 'rand', IF_VERBOSE = True, kappa = 0):
    N = NE
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    if IF_VERBOSE:
    	print baseFldr
    filename = 'MI1_inst_theta%.6f_tr%s.txt'%(phi, trNo)
    return np.loadtxt(baseFldr + filename, delimiter = ';')

def PlotM1vsT(p = 0, gamma = 0, phi = 0, trNo = 0, mExt = 0.075, mExtOne = 0.075, K = 1000, NE = 10000, T = 1000, nPop = 2, rewireType = 'cntrl', IF_VERBOSE = True, kappa = 0):
    out = LoadM1vsT(p, gamma, phi, trNo, mExt, mExtOne, K, NE, T, nPop, rewireType, IF_VERBOSE, kappa)
    _, nColumns = out.shape
    if nColumns == 3:
        m1 = out[:, 0]
        m1Phase = out[:, 1]
        phi_ext = out[:, 2]
        tAxis = np.linspace(0, 1, m1.size)
        plt.plot(tAxis, m1)
        plt.xlabel('Time (a.u)')
        plt.ylabel(r'$m_E^{(1)}$')
        plt.ylim(0, .25)
        plt.vlines(.3, *plt.ylim(), color = 'k')
        plt.vlines(.7, *plt.ylim(), color = 'k')        
        filename = './figs/twopop/' + 'stimulus_tracking_m1_p%sgmma%s'%(int(p*10), int(gamma *100))
	paperSize = [2.5, 2]
        axPosition = [0.22, 0.2, .7, .7]
        print 'printing figure'
        ProcessFigure(plt.gcf(), filename, 1, paperSize = paperSize, axPosition=axPosition, titleSize=10, nDecimalsX=1, nDecimalsY=1, figFormat='pdf', labelFontsize = 8, tickFontsize = 6)
        plt.figure()
        plt.plot(tAxis, phi_ext*180/np.pi, '--', label = r'stimulus')    
        plt.plot(tAxis, m1Phase*180/np.pi, alpha = 0.5, label = r'$\angle m_E(\phi)$')
        plt.ylim(0, 180)
        plt.xlabel('Time (a.u)')
        plt.ylabel(r'Phase (deg)')
        plt.legend(loc = 2, frameon = False, numpoints = 1, ncol = 1, prop = {'size': 8})            
        filename = './figs/twopop/' + 'stimulus_tracking_phase_p%sgmma%s'%(int(p*10), int(gamma *100))
	paperSize = [2.5, 2]
        axPosition = [0.22, 0.2, .7, .7]
        print 'printing figure'
        ProcessFigure(plt.gcf(), filename, 1, paperSize = paperSize, axPosition=axPosition, titleSize=10, nDecimalsX=1, nDecimalsY=1, figFormat='pdf', labelFontsize = 8, tickFontsize = 6)        
        ipdb.set_trace()
    else:
        print '-' * 25
        print 'no simulus change'
        print '-' * 25        

def M1Component(x):
    out = np.nan
    if len(x) > 0:
	dPhi = np.pi / len(x)
	out = 2.0 * np.absolute(np.dot(x, np.exp(-2.0j * np.arange(len(x)) * dPhi))) / len(x)
    return out

def M2Component(x):
    out = np.nan
    if len(x) > 0:
	dPhi = np.pi / len(x)
	out = 2.0 * np.absolute(np.dot(x, np.exp(-4.0j * np.arange(len(x)) * dPhi))) / len(x)
    return out

def PopM2Component(tc):
    nNeurons, nPhis = tc.shape
    outM2 = np.zeros((nNeurons, ))
    for i in range(nNeurons):
        outM2[i] = M2Component(tc[i, :])
    return outM2

def PopM1Component(tc):
    nNeurons, nPhis = tc.shape
    outM1 = np.zeros((nNeurons, ))
    for i in range(nNeurons):
        outM1[i] = M1Component(tc[i, :])
    return outM1    

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
	    mr = LoadFr(p, gamma, phi, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, IF_VERBOSE = IF_VERBOSE, kappa = kappa)
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

def KappaVsM1(kappaList, nTrials = 10, p=0, gamma=0, nPhis = 8, mExt=0.075, mExtOne=0.075, rewireType='rand', N=10000, K=1000, nPop=2, T=1000, IF_PO_SORTED = False, sortedIdx = [], minRate=0):
    m1E = np.zeros((nTrials + 2, len(kappaList)))
    m1E[0, :] = np.nan
    m0E = np.zeros((nTrials + 2, len(kappaList)))
    m0E[0, :] = np.nan    
    for trNo in range(1, nTrials + 2): #$ trNo 0 is always the CONTROL 
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

def KappaVsM1overM0(kappaList, nTrials = 10, p=0, gamma=0, nPhis = 8, mExt=0.075, mExtOne=0.075, rewireType='rand', N=10000, K=1000, nPop=2, T=1000, IF_PO_SORTED = False, sortedIdx = [], minRate = 0):
    m1E = np.zeros((nTrials + 2, len(kappaList)))
    m1E[0, :] = np.nan
    m0E = np.zeros((nTrials + 2, len(kappaList)))
    m0E[0, :] = np.nan    
    for trNo in range(1, nTrials + 2): #$ trNo 0 is always the CONTROL 
	print ''
	print '--' * 27
	print 'tr#: ', trNo
	print '--' * 27
        # m1E[trNo, :] = KappaVsM1AtTr(kappaList, p, gamma, nPhis, mExt, mExtOne, rewireType, N, K, nPop, T, trNo, IF_PO_SORTED, sortedIdx)
        m1E[trNo, :], m0E[trNo, :] = KappaVsM1AtTr(kappaList, p, gamma, nPhis, mExt, mExtOne, rewireType, N, K, nPop, T, trNo, IF_PO_SORTED, sortedIdx, minRate)
    # m1EvsKappa = np.nanmean(m1E, 0)
    # m0EvsKappa = np.nanmean(m0E, 0)    
    nValidTrials = np.sum(~np.isnan(m1E.mean(1)))
    m1overM0 = m1E / m0E
    if nValidTrials >= 2:
	m1EvsKappaSEM = np.nanstd(m1overM0, 0) / np.sqrt(float(nValidTrials))
    else:
	m1EvsKappaSEM = np.nan
    (_, caps, _) = plt.errorbar(kappaList, np.nanmean(m1overM0, 0), fmt = 'o-', markersize = 3, yerr = m1EvsKappaSEM, lw = 0.8, elinewidth=0.8, label = r'$N = %s, K = %s$'%(N, K))
    for cap in caps:
        cap.set_markeredgewidth(0.8)
    plt.xlim(min(kappaList) - 1, max(kappaList) + 1)
    plt.xlabel(r'$\kappa$', fontsize = 20)
    plt.ylabel(r'$\frac{m_E^{(1)}}{m_E^{(0)}}$', fontsize = 20)
    plt.show()
    ipdb.set_trace()

    
    # return m1EvsKappa, m1EvsKappaSEM

        

def CompareMeanOSIvsKappa(kappaList, mExtOne, p = 0, gamma = 0, nPhis = 8, mExt = 0.075,  trNo = 1, rewireType = 'rand', N = 10000, K = 1000, nPop = 2, T = 1000, IF_NEW_FIG = True, clrCntr = 0, filename = '', IF_LEGEND = True, legendTxt = '', color = '', neuronType = 'E'):
    # colors = [plt.cm.Dark2(i) for i in np.linspace(0, 1, 1 + clrCntr + len(pList) * len(gList) * len(mExtOneList) * len(trList), endpoint = False)]
    # colors = [plt.cm.Dark2(i) for i in np.linspace(0, 1, clrCntr + len(kappaList) * len(gList) * len(mExtOneList) * len(KList), endpoint = False)]
    colors = [plt.cm.Dark2(i) for i in np.linspace(0, 1, clrCntr + len(kappaList), endpoint = False)]    
    meanOSI = []
    meanOSI_I = []
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    validTrList = []
    legendTxtTag = legendTxt
    if IF_NEW_FIG:
	plt.figure()
    validTrList = []
    meanOSI = []
    meanOSI_I = []
    counter = 0
    validKappa = []
    for kk, kappa in enumerate(kappaList):
	# legendTxt = ', K=%s'%(K) + legendTxtTag
	legendTxt = ', %s'%(neuronType) + legendTxtTag
	try:
	    print trNo
	    tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
	    osi = OSIOfPop(tc[:N, :], phis)
	    osiI = OSIOfPop(tc[N:, :], phis)
	    meanOSI.append(np.nanmean(osi))
	    meanOSI_I.append(np.nanmean(osiI))
	    if ~np.isnan(meanOSI[counter]):
		validTrList.append(trNo)
	    clrCntr += 1
	    counter += 1
	    validKappa.append(kappa)
	except IOError:
	    print "p = ", p, " gamma = ", gamma, " trial# ", trNo, " file not found"
    meanOSI = np.array(meanOSI)
    meanOSI_I = np.array(meanOSI_I)
    kappaList = np.asarray(kappaList, dtype = float)
    validIdx = ~np.isnan(meanOSI)
    validKappa = kappaList[validIdx]
    pcolor = color
    if color == '':
	pcolor = colors[kk]
    if neuronType == 'E':
	plt.plot(validKappa, meanOSI[validIdx], 'o-', color = pcolor, label = r'$m_0^{(1)}=%s$'%(mExtOne) + legendTxt, markeredgecolor = pcolor)
    else:
	plt.plot(validKappa, meanOSI_I[validIdx], 'o-', color = pcolor, label = r'$m_0^{(1)}=%s$'%(mExtOne) + legendTxt, markeredgecolor = pcolor)
    print '--'*26
    osiLast = meanOSI[-1]
    osilastCnt = -1
    osiLastI = meanOSI_I[-1]
    while np.isnan(osiLast):
	osilastCnt -= 1
	print osilastCnt
	osiLast = meanOSI[osilastCnt]	
	osiLastI = meanOSI_I[osilastCnt]
    print 'pc change in mean OSI = ', 100 * (osiLast - meanOSI[0]) / meanOSI[0]
    print '--'*26
    plt.gca().set_position([0.15, 0.15, .65, .65])
    plt.ylim(0, .5)
    xmax = max(validKappa) + 1
    xmin = min(validKappa) - 1
    plt.xlim(xmin, xmax)
    plt.xlim(0, 12)    
    plt.ylim(0, 1)
    plt.xlabel(r'$\kappa$')
    plt.ylabel(r'$\langle OSI \rangle$')
    plt.gca().set_position([0.25, 0.25, .65, .65])
    filename = filename + "p%sg%sk%s_K%s_"%(p, gamma, kappa, K) + neuronType
    filepath = "~/binary/figs/rewiring/compareOSI_mean_"+ rewireType + '_' +  filename
    print filepath
    paperSize = [4, 3]
    plt.legend(loc = 2, frameon = False, numpoints = 1, prop = {'size': 8})
    print 'saving as: ', filepath
    ProcessFigure(plt.gcf(), filename, IF_SAVE = 1, IF_XTICK_INT = False, figFormat = 'pdf')
    plt.show()
    print meanOSI

def CumulativeDistr(x):
    x.sort()
    # F2 = np.array(range(N))/float(N)    
    plt.step(np.concatenate([x, x[[-1]]]), np.arange(x.size+1) / float(x.size + 1))
    # plt.plot(x, F2)

def CompareMeanOSIvsKappaRewire(kappaList, nTrials = 10, mExtOne=0.075, p = 0, gamma = 0, nPhis = 8, mExt = 0.075, rewireType = 'rand', N = 10000, K = 1000, nPop = 2, T = 1000, IF_NEW_FIG = True, clrCntr = 0, filename = '', IF_LEGEND = True, legendTxt = '', color = '', neuronType = 'E', markerType = 'o-'):

    osiE = np.empty((len(kappaList), nTrials))
    osiI = np.empty((len(kappaList), nTrials))
    osiE[:] = np.nan
    osiI[:] = np.nan
    osiESEM = np.empty((len(kappaList), ))
    osiISEM = np.empty((len(kappaList), ))
    osiESEM[:] = np.nan
    osiISEM[:] = np.nan

    # trList = range(3000, 3000 + nTrials)
    trList = range(2000, 2000 + nTrials)    
    trListOld = trList
    trList = range(100, 100 + nTrials)    
    mE = np.empty((len(kappaList), nTrials))
    mI = np.empty((len(kappaList), nTrials))
    mESEM = np.empty((len(kappaList), ))
    mISEM = np.empty((len(kappaList), ))
    mESEM[:] = np.nan; mISEM[:] = np.nan
    for idx, kappa in enumerate(kappaList):
        nValidTrials = 0
	# if kappa == 0 and K == 2000:
	#     trList = [0]
	if kappa == 8 and K == 1000:
	    trList = trListOld
        for trIdx, trNo in enumerate(trList):
            print 'trIdx', trIdx
            # ipdb.set_trace()
            tc, IF_FILE_LOADED = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, kappa=kappa, trNo = trNo, IF_SUCCESS = True, K=K)
            mE[idx, trIdx] = np.nanmean(tc[:N]); mI[idx, trIdx] = np.nanmean(tc[N:])
            tmp = OSIOfPop(tc, np.pi)
            osiE[idx, trIdx] = np.nanmean(tmp[:N])
            osiI[idx, trIdx] = np.nanmean(tmp[N:])
            nValidTrials += IF_FILE_LOADED
        osiESEM[idx] = np.nanstd(osiE[idx, :]) / np.sqrt(nValidTrials)
        osiISEM[idx] = np.nanstd(osiI[idx, :]) / np.sqrt(nValidTrials)
        mESEM[idx] = np.nanstd(mE[idx, :]) / np.sqrt(nValidTrials)
        mISEM[idx] = np.nanstd(mI[idx, :]) / np.sqrt(nValidTrials)        

    # plt.plot(osiE.mean(1), 'ko')
    # plt.plot(osiI.mean(1), 'ro')    
    # ipdb.set_trace()
    # osiESEM = np.nanstd(osiE, 1) / np.sqrt(nValidTrials)
    # osiISEM = np.nanstd(osiI, 1) / np.sqrt(nValidTrials)    

    plt.figure(0)
    (_, caps, _) = plt.errorbar(range(len(kappaList)), np.nanmean(osiE, 1), color = 'k', fmt = markerType, markersize = 4, yerr = osiESEM, lw = 0.8, elinewidth=0.8, label = 'E, ' + legendTxt, markeredgecolor = 'k')
    for cap in caps:
        cap.set_markeredgewidth(0.8)
    (_, caps, _) = plt.errorbar(range(len(kappaList)), np.nanmean(osiI, 1), color = 'r', fmt = markerType, markersize = 4, yerr = osiISEM, lw = 0.8, elinewidth=0.8, label = 'I, ' + legendTxt, markeredgecolor = 'r')
    for cap in caps:
        cap.set_markeredgewidth(0.8)
 
    plt.xlabel(r'$\kappa$')
    plt.ylabel(r'$\langle OSI \rangle$')
    plt.legend(loc = 2, frameon = False, numpoints = 1, prop = {'size': 6})
    plt.xlim(-.5, 1.5)
    plt.ylim(0, .5)

    filename = './figs/PNAS/twopop/' + 'meanOSI_vs_kappa_rewired_vsK_tr%s'%(trList[0])
    paperSize = [2.5, 2]
    axPosition=[.24, .24, .75, .7]
    ProcessFigure(plt.gcf(), filename, 1, paperSize = paperSize, axPosition = axPosition, titleSize=10, nDecimalsX=1, nDecimalsY=2, figFormat='svg', labelFontsize = 8, tickFontsize = 6)

    plt.xticks([0, 1])
    plt.gca().set_xticklabels(kappaList)
    Print2Pdf(plt.gcf(),  filename,  paperSize, figFormat='svg', labelFontsize = 8, tickFontsize=6, titleSize = 10.0)
    plt.show()


    plt.figure(1)
    (_, caps, _) = plt.errorbar(range(len(kappaList)), np.nanmean(mE, 1), color = 'k', fmt = markerType, markersize = 4, yerr = mESEM, lw = 0.8, elinewidth=0.8, label = 'E, ' + legendTxt, markeredgecolor = 'k')
    for cap in caps:
        cap.set_markeredgewidth(0.8)
    (_, caps, _) = plt.errorbar(range(len(kappaList)), np.nanmean(mI, 1), color = 'r', fmt = markerType, markersize = 4, yerr = mISEM, lw = 0.8, elinewidth=0.8, label = 'I, ' + legendTxt, markeredgecolor = 'r')
    for cap in caps:
        cap.set_markeredgewidth(0.8)
    plt.xlabel(r'$\kappa$')
    plt.ylabel(r'$\langle m_A \rangle$')
    plt.ylim(0, 1)
    plt.xlim(-.5, 1.5)
    plt.legend(loc = 2, frameon = False, numpoints = 1, prop = {'size': 6})        
    filename = './figs/PNAS/twopop/' + 'meanRates_vs_kappa_rewired_vsK_trNo%s'%(trList[0])
    paperSize = [2.5, 2]
    axPosition=[.24, .24, .75, .7]
    ProcessFigure(plt.gcf(), filename, 1, paperSize = paperSize, axPosition = axPosition, titleSize=10, nDecimalsX=1, nDecimalsY=2, figFormat='svg', labelFontsize = 8, tickFontsize = 6)
    plt.xticks([0, 1])
    plt.gca().set_xticklabels(kappaList)
    Print2Pdf(plt.gcf(),  filename,  paperSize, figFormat='svg', labelFontsize = 8, tickFontsize=6, titleSize = 10.0)
    plt.show()


    print 'K = ', K
    print 'osi: E', np.nanmean(osiE, 1)
    print 'ratio: ', np.nanmean(osiE, 1)[1] /  np.nanmean(osiE, 1)[0]
    print 'osi: I', np.nanmean(osiI, 1)
    print 'ratio: ', np.nanmean(osiI, 1)[1] /  np.nanmean(osiI, 1)[0]

    print 'm: E', np.nanmean(mE, 1)
    print 'ratio: ', np.nanmean(mE, 1)[1] /  np.nanmean(mE, 1)[0]
    print 'm: I', np.nanmean(mI, 1)
    print 'ratio: ', np.nanmean(mI, 1)[1] /  np.nanmean(mI, 1)[0]    
    
    # ipdb.set_trace()    
    
    
    
    
def ComputeAllCCC(kappa, trNo, p=0, gamma=0, nPhis=8, mExt=0.075, mExtOne=0.075, rewireType='rand', N=10000, K=1000, nPop=2, T=1000, NFF=10000, JE0 = 2.0, JI0 = 1.0, IF_COMPUTE = False):
    out = np.nan
    if IF_COMPUTE:
	try:
	    tcIn = GetInputTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa, 'net')
	    tcOut = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
	    theta = np.linspace(0, 180, nPhis, endpoint = False)
	    poIn = POofPopulation(tcIn, theta, IF_IN_RANGE = True) * np.pi / 180
	    poOut = POofPopulation(tcOut, theta, IF_IN_RANGE = True) * np.pi / 180
	    ##################################################
	    uFF = rw.ComputeFFInput(nPhis, p, gamma, kappa, mExt, mExtOne, trNo, N, K, nPop, NFF, JE0, JI0, 0.2, rewireType, T)
	    poFF = POofPopulation(uFF, theta, IF_IN_RANGE = True) * np.pi / 180
	    tcOut = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
	    poOut = POofPopulation(tcOut, theta, IF_IN_RANGE = True) * np.pi / 180
	    print 'computing in-FF ccc... ', trNo
	    sys.stdout.flush()
	    out = CircularCorrCoeff(poFF[:N], poOut[:N])
	    np.save('./data/NET_IN_FF_CCC_kappa%s_tr%s_K%s'%(int(kappa*10), int(trNo), int(K)), out)                
	    print 'done', trNo, 'CCC=', out
	    ##################################################
	    print 'computing in-out ccc... ', trNo
	    sys.stdout.flush()
	    outInOut = CircularCorrCoeff(poIn[:N], poOut[:N])
	    np.save('./data/NET_IN_OUT_CCC_kappa%s_tr%s_K%s'%(int(kappa*10), int(trNo), K), outInOut)        
	    print 'done', trNo
	    ##################################################
	    uRec = tcIn - uFF
	    poRec = POofPopulation(uRec, theta, IF_IN_RANGE = True) * np.pi / 180
	    print 'computing rec-FF ccc... ', trNo
	    sys.stdout.flush()
	    out = CircularCorrCoeff(poRec[:N], poOut[:N])
	    np.save('./data/REC_IN_FF_CCC_kappa%s_tr%s_K%s'%(int(kappa*10), int(trNo), K), out)                        
	    print 'done', trNo, 'CCC=', out
	    ##################################################
	except IOError:
	    print ''
	return out
    else:
	if K == 1000:
	    kappa = 0; trNo = 0
	    cccNET_FF_Cntr = np.load('./data/NET_IN_FF_CCC_kappa%s_tr%s_'%(int(kappa*10), int(trNo)) + '.npy')
	    cccNET_OUT_Cntr = np.load('./data/NET_IN_OUT_CCC_kappa%s_tr%s_'%(int(kappa*10), int(trNo)) + '.npy')
	    cccREC_FF_Cntr = np.load('./data/REC_IN_FF_CCC_kappa%s_tr%s_'%(int(kappa*10), int(trNo)) + '.npy')

	    kappa = 8; trNo = 100
	    cccNET_FF_k8 = np.load('./data/NET_IN_FF_CCC_kappa%s_tr%s_'%(int(kappa*10), int(trNo)) + '.npy')
	    cccNET_OUT_k8 = np.load('./data/NET_IN_OUT_CCC_kappa%s_tr%s_'%(int(kappa*10), int(trNo)) + '.npy')
	    cccREC_FF_k8 = np.load('./data/REC_IN_FF_CCC_kappa%s_tr%s_'%(int(kappa*10), int(trNo)) + '.npy')
	else:
	    ipdb.set_trace()
	    kappa = 0; trNo = 0
	    cccNET_FF_Cntr = np.load('./data/NET_IN_FF_CCC_kappa%s_tr%s_K%s'%(int(kappa*10), int(trNo), K) + '.npy')
	    cccNET_OUT_Cntr = np.load('./data/NET_IN_OUT_CCC_kappa%s_tr%s_K%s'%(int(kappa*10), int(trNo), K) + '.npy')
	    cccREC_FF_Cntr = np.load('./data/REC_IN_FF_CCC_kappa%s_tr%s_K%s'%(int(kappa*10), int(trNo), K) + '.npy')

	    kappa = 8; trNo = 100
	    cccNET_FF_k8 = np.load('./data/NET_IN_FF_CCC_kappa%s_tr%s_K%s'%(int(kappa*10), int(trNo), K) + '.npy')
	    cccNET_OUT_k8 = np.load('./data/NET_IN_OUT_CCC_kappa%s_tr%s_K%s'%(int(kappa*10), int(trNo), K) + '.npy')
	    cccREC_FF_k8 = np.load('./data/REC_IN_FF_CCC_kappa%s_tr%s_K%s'%(int(kappa*10), int(trNo), K) + '.npy')
	    
	    

        # ipdb.set_trace() 
        xax = [1, 2]
	plt.plot(xax, [cccNET_OUT_Cntr, cccNET_OUT_k8], 'go-', label = 'out-net')
	plt.plot(xax, [cccNET_FF_Cntr, cccNET_FF_k8], 'ko-', label = 'out-FF')
	plt.plot(xax, [cccREC_FF_Cntr, cccREC_FF_k8], 'co-', label = 'out-rec')

	plt.xlabel(r'$\kappa$')
	plt.ylabel('Circ Corr Coeff')
	plt.ylim(0, 1)
	plt.xlim(0.5, 2.5)
	plt.xticks([1, 2])
	plt.gca().set_xticklabels([0, 8])
	plt.show()

def PlotTuningCurves(tc, filetag):
    _, nThetas = tc.shape
    thetas = np.arange(-90, 91, 180./8)

    neuronsIdx = [1798, 7753, 3233, 6218]
    print neuronsIdx
    ###### EXCITATORY 
    for idx in neuronsIdx:
        tmp = tc[idx, :]
        # tmp = tmp / np.max(tmp)
        tmp = np.roll(tmp, -1 * np.argmax(tmp))
        tmp = np.roll(tmp, 4)        
        plt.plot(thetas, np.concatenate((tmp, [tmp[0]])), 'o-k', lw = .5, markersize = 0.95)
        filename = './figs/PNAS/tc/tuning_curves_%s_idx_%s'%(filetag, idx)
        paperSize = [2.5/2, 2/2.0]
        axPosition=[.25, .25, .65, .65]
        FixAxisLimits(plt.gcf())
        plt.gca().set_xticks([-90, 0, 90])
        ymin, ymax = plt.ylim()
        plt.ylim(0, ymax)
        plt.gca().set_yticks([0, ymax * 0.5, ymax])
        plt.gca().set_yticklabels(['0', '%.4s'%(ymax * 0.5), '%.4s'%(ymax)])
        plt.gca().set_xticks([0, 90, 180])
        plt.gca().set_xticklabels(['0', '90', '180'])        
        plt.draw()
        ProcessFigure(plt.gcf(), filename, 1, paperSize = paperSize, axPosition = axPosition, titleSize=10, nDecimalsX=1, nDecimalsY=1, figFormat='svg', labelFontsize = 6, tickFontsize = 6)
        plt.clf()
    ##### INHIBITORY
    NE = 10000
    neuronsIdx = np.random.choice(10000, 20) + NE
    neuronsIdx = [13023, 14079, 18808, 18853]
    print neuronsIdx
    for idx in neuronsIdx:
        tmp = tc[idx, :]
        # tmp = tmp / np.max(tmp)
        tmp = np.roll(tmp, -1 * np.argmax(tmp))
        tmp = np.roll(tmp, 4)        
        plt.plot(thetas, np.concatenate((tmp, [tmp[0]])), 'o-r', lw = .5, markersize = 0.95)
        filename = './PUB_FIGS/tuning_curves_%s_I_idx_%s'%(filetag, idx)
        paperSize = [2.5/2, 2/2.0]
        axPosition=[.25, .25, .65, .65]
        FixAxisLimits(plt.gcf())
        plt.gca().set_xticks([-90, 0, 90])
        plt.draw()
        ProcessFigure(plt.gcf(), filename, 1, paperSize = paperSize, axPosition = axPosition, titleSize=10, nDecimalsX=1, nDecimalsY=1, figFormat='svg', labelFontsize = 6, tickFontsize = 6)
        plt.clf()
    

def ComapreTuningCurves(kappaList, nTrials = 1, mExtOne=0.075, p = 0, gamma = 0, nPhis = 8, mExt = 0.075, rewireType = 'rand', N = 10000, K = 1000, nPop = 2, T = 1000, IF_NEW_FIG = True, clrCntr = 0, filename = '', IF_LEGEND = True, legendTxt = '', color = '', neuronType = 'E', markerType = 'o-'):
    tclist = []

    trListOld = range(2000, 2000 + nTrials)
    # trList = range(100, 100 + nTrials)


    
    trList = range(100, 100 + nTrials)
    
    for idx, kappa in enumerate(kappaList):
        nValidTrials = 0
	# if kappa == 0 and K == 2000:
	#     trList = [0]
	if kappa == 8 and K == 1000:
	    trList = trListOld
        for trIdx, trNo in enumerate(trList):
            print 'trIdx', trIdx
            # ipdb.set_trace()
            tc, IF_FILE_LOADED = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, kappa=kappa, trNo = trNo, IF_SUCCESS = True, K=K)
        if IF_FILE_LOADED:
            tclist.append(tc)
    assert len(tclist) == len(kappaList)
    PlotTuningCurvesKappaNonZero(tclist, kappaList, neuronType = neuronType)

        
def PlotTuningCurvesKappaNonZero(tcList, kappaList, filetag='', neuronType='E', roll=False):
    _, nThetas = tcList[0].shape
    # thetas = np.arange(-90, 91, 180./8)
    thetas = np.arange(0., 181., 180./8)    
    neuronsIdx = [1798, 7753, 3233, 6218]
    if neuronType == 'E':
        neuronsIdx = np.random.choice(10000, 100)
    else:
        neuronsIdx = np.random.choice(10000, 100) + 10000

    print neuronsIdx
    filetag = neuronType
    ###### EXCITATORY
    colors = ['k', 'g', 'r']
    for idx in neuronsIdx:    
        for kappaIdx, kappa in enumerate(kappaList):
            tc = tcList[kappaIdx]
            tmp = tc[idx, :]
            if roll:
                tmp = np.roll(tmp, -1 * np.argmax(tmp))
                tmp = np.roll(tmp, 4)        
            plt.plot(thetas, np.concatenate((tmp, [tmp[0]])), 'o-', color = colors[kappaIdx], lw = .5, markersize = 0.95, markerfacecolor = colors[kappaIdx], markeredgecolor=colors[kappaIdx])
        filename = './figs/PNAS/tc/tuning_curves_%s_idx_%s'%(filetag, idx)
        paperSize = [2.5/2, 2/2.0]
        axPosition=[.25, .25, .65, .65]
        FixAxisLimits(plt.gcf())

        ymin, ymax = plt.ylim()
        plt.ylim(0, ymax)
        plt.gca().set_yticks([0, ymax * 0.5, ymax])
        plt.gca().set_yticklabels(['0', '%.4s'%(ymax * 0.5), '%.4s'%(ymax)])
        plt.gca().set_xticks([0, 90, 180])
        plt.gca().set_xticklabels(['0', '90', '180'])        
        plt.draw()
        ProcessFigure(plt.gcf(), filename, 1, paperSize = paperSize, axPosition = axPosition, titleSize=10, nDecimalsX=1, nDecimalsY=1, figFormat='svg', labelFontsize = 6, tickFontsize = 6)
        plt.clf()


def PlotAvgOSI(kappaList, NE = 10000):
    
    tckappa = []
    osiEList = []
    osiIList = []
    validIdx = []
    
    for kappaIdx, kappa in enumerate(kappaList):
        try:
            if kappa == 0:
                tckappa = GetTuningCurves(0, 0, 8, 0.075, 0.075, 'rand', kappa=0)
            else:
                tckappa = GetTuningCurves(0, 0, 8, 0.075, 0.075, 'rand', kappa=kappa, trNo=2)
            osiEList.append(np.nanmean(OSIOfPop(tckappa[:NE], 0)))
            osiIList.append(np.nanmean(OSIOfPop(tckappa[NE:], 0)))
            validIdx.append(kappaIdx)
        except IOError:
            pass
    print kappaList[validIdx]
    plt.plot(kappaList[validIdx], osiEList)
    plt.plot(kappaList[validIdx], osiIList)    
    
                          
    

def PrintTuningBook(tcCntr, tcNew, nNeurons, fname, neuronType='E', NE = 10000, NI = 10000, nPhis = 8, color = 'k', kappaNew = 8):
    doc = Document(fname)
    doc.packages.append(Package('geometry', options=['left=2cm', 'right=2cm']))
    nFigsPerPage = 12
    nPages = int(np.ceil(nNeurons / float(nFigsPerPage)))
    plt.figure()
    plt.ioff()
    theta = np.linspace(0, 180, nPhis, endpoint = False)
    doc = Document(fname)
    doc.packages.append(Package('geometry', options=['left=2cm', 'right=2cm']))
    width = r'.3\linewidth'
    osiCntr = OSIOfPop(tcCntr, np.pi)
    osiNew = OSIOfPop(tcNew, np.pi)    
    if neuronType == 'E':
	rndNeurons = np.random.randint(0, NE, nNeurons)
    else:
	rndNeurons = np.random.randint(NE, NE + NI, nNeurons)	
    for kk in range(nPages):
	rndNeuronsx = rndNeurons[kk * nFigsPerPage: (kk + 1) * nFigsPerPage]
	print kk * nFigsPerPage, (kk + 1) * nFigsPerPage
	with doc.create(Figure(position='htbp')) as plot:
	    for idx, i in enumerate(rndNeuronsx):
                labelTxtCntrl = r'$\kappa = 0, \, \mathrm{OSI} = %.4s $'%(osiCntr[i])
                labelTxtNew = r'$\kappa = %s, \, \mathrm{OSI} = %.4s $'%(kappaNew, osiNew[i])
		plt.plot(np.concatenate((theta, [180])), np.concatenate((tcCntr[i, :], [tcCntr[i, 0]])), 'ko-', label = labelTxtCntrl)
		plt.plot(np.concatenate((theta, [180])), np.concatenate((tcNew[i, :], [tcNew[i, 0]])), 'o-', color = color, label = labelTxtNew, markeredgecolor = color)
		mainAxis = plt.gca()
		mainAxis.set_title('neuron#%s'%(i), fontsize = 16)
                plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 12})		
		with doc.create(SubFigure(position='b', width=NoEscape(width))) as figure:
		    figure.add_plot(width=NoEscape(r'\linewidth'), dpi = 300) #*args, **kwargs)
		plt.clf()
    doc.generate_pdf(clean_tex=False)

def PrintInputTuningBook(a, b, c, d, e, nNeurons, fname, neuronType='E', NE = 10000, NI = 10000, nPhis = 8, color = 'k', kappaNew = 8):
    doc = Document(fname)
    doc.packages.append(Package('geometry', options=['left=2cm', 'right=2cm']))
    nFigsPerPage = 12
    nPages = int(np.ceil(nNeurons / float(nFigsPerPage)))
    plt.figure()
    plt.ioff()
    theta = np.linspace(0, 180, nPhis, endpoint = False)
    doc = Document(fname)
    doc.packages.append(Package('geometry', options=['left=2cm', 'right=2cm']))
    width = r'.3\linewidth'
    osiCntr = OSIOfPop(a, np.pi)
    osiNew = OSIOfPop(b, np.pi)
    d = a[:NE, :] + b[:NE, :] + c[:NE, :]
    if neuronType == 'E':
	rndNeurons = np.random.randint(0, NE, nNeurons)
    else:
	rndNeurons = np.random.randint(NE, NE + NI, nNeurons)	
    for kk in range(nPages):
	rndNeuronsx = rndNeurons[kk * nFigsPerPage: (kk + 1) * nFigsPerPage]
	print kk * nFigsPerPage, (kk + 1) * nFigsPerPage
	with doc.create(Figure(position='htbp')) as plot:
	    for idx, i in enumerate(rndNeuronsx):
                labelTxtCntrl = r'$\kappa = 0, \, \mathrm{OSI} = %.4s $'%(osiCntr[i])
                labelTxtNew = r'$\kappa = %s, \, \mathrm{OSI} = %.4s $'%(kappaNew, osiNew[i])
		labelTxtA = r'$u_{EE}$'
		labelTxtB = r'$abs(u_{EI})$'
		labelTxtC = r'$u_{FF}$'
		labelTxtD = r'$u_{net}$'
		labelTxtE = r'$m(\phi)$'
		# plt.plot(np.concatenate((theta, [180])), np.concatenate((a[i, :], [a[i, 0]])) / np.max(a[i, :]), 'ko-', label = labelTxtA)
		# plt.plot(np.concatenate((theta, [180])), np.concatenate((b[i, :], [b[i, 0]])) / np.max(np.abs(b[i, :])), 'bo-', label = labelTxtB, markeredgecolor = 'b')
		plt.plot(np.concatenate((theta, [180])), np.concatenate((b[i, :], [b[i, 0]])), 'bo-', label = labelTxtB, markeredgecolor = 'b')		
		# plt.plot(np.concatenate((theta, [180])), np.concatenate((c[i, :], [c[i, 0]])) / np.abs(np.max(c[i, :])), 'go-', label = labelTxtC, markeredgecolor = 'g')
		# plt.plot(np.concatenate((theta, [180])), np.abs(np.concatenate((d[i, :], [d[i, 0]]))) / np.max(np.abs(d[i, :])), 'co-', label = labelTxtD, markeredgecolor = 'c')
		# plt.plot(np.concatenate((theta, [180])), np.concatenate((e[i, :], [e[i, 0]])) / np.max(e[i, :]), 'ro-', label = labelTxtE, markeredgecolor = 'r')
		mainAxis = plt.gca()
		mainAxis.set_title('neuron#%s'%(i), fontsize = 16)
		# plt.ylim(, 1.05)
                plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 12})		
		with doc.create(SubFigure(position='b', width=NoEscape(width))) as figure:
		    figure.add_plot(width=NoEscape(r'\linewidth'), dpi = 300) #*args, **kwargs)
		plt.clf()
    doc.generate_pdf(clean_tex=False)
    
def PlotMeanTc(tc, kappa, J, nPhis = 8, NE = 10000, labelTxt='', pcolor='k'):
    # J is the rewired strenghtened prefactor
    # thetas = np.linspace(0, 180, nPhis, endpoint = 1)
    thetas = np.arange(-90, 90, 180.0/8)
    prefferedOri = np.argmax(tc, 1)
    tcmax = np.max(np.abs(tc), 1)
    # ipdb.set_trace()
    tcmax.shape = NE, 1
    tcmax = np.tile(tcmax, (1, nPhis))
    tc = tc / tcmax
    cvMat = np.empty((NE, len(thetas)))
    for kNeuron in np.arange(NE):
	cvMat[kNeuron, :] = np.roll(tc[kNeuron, :], -1 * prefferedOri[kNeuron])
    plt.ion()
    tmpE = cvMat #[plotId[plotId < NE], :]
    meanE = np.nanmean(tmpE, 0)
    meanE = np.roll(meanE, 4)
    osi = OSI(meanE, np.arange(0, 180, 22.5))
    print 'osi = ', osi
    thetas = np.arange(-90, 91, 22.5)
    plt.plot(thetas, np.concatenate((meanE, [meanE[0]])), 'o-', label=labelTxt + ' osi: %.5s'%(osi), color = pcolor, markersize = 2, markeredgecolor = pcolor)
    # plt.ylim(0.8, 1)
    plt.gca().set_xticks([-90, 0, 90])

    filename = './figs/twopop/rewire/' + 'pop_mean_tuning_kappa%s'%(kappa)
    paperSize = [2.5, 2]
    axPosition=[.22, .22, .65, .65]
    plt.ylim([0, 1])
    plt.gca().set_yticks([0, 0.5, 1])
    ProcessFigure(plt.gcf(), filename, 1, paperSize = paperSize, axPosition = axPosition, titleSize=10, nDecimalsX=1, nDecimalsY=1, figFormat='svg', labelFontsize = 10, tickFontsize = 8)
    
    


def PlotMeanUEUITuning(ue, ui, kappa, J, IF_LEGEND = 0):
    PlotMeanTc(ue, kappa, J, labelTxt='EE', pcolor='k')
    PlotMeanTc(ui, kappa, J, labelTxt='EI', pcolor='r')
    plt.title(r'$\kappa = %s, \, \omega = %s$'%(kappa, J))
    plt.xlabel('PO (deg)')
    plt.ylabel('Input')
    if IF_LEGEND:
        plt.legend(loc = 8, frameon = False, numpoints = 1, prop = {'size': 6})
    filename = './figs/twopop/rewire/' + 'mean_rec_input_kappa%s_%sJ'%(kappa, J)
    paperSize = [2.5, 2]
    axPosition=[.22, .22, .65, .65]
    ProcessFigure(plt.gcf(), filename, 1, paperSize = paperSize, axPosition = axPosition, titleSize=10, nDecimalsX=1, nDecimalsY=1, figFormat='pdf', labelFontsize = 10, tickFontsize = 8)
    

def PlotUePO(idxvec, sparsevec, npost, tc, JEE = 1, JEI = -1.5, JIE = 1, JII = -1, NE = 10000, NI = 10000, nPhis = 8):
    fr = tc.mean(1) # mean over angles
    po = POofPopulation(tc, IF_IN_RANGE = 1)
    nNeurons = idxvec.size
    ue = np.zeros((nNeurons, nPhis))
    for i in range(NE):
	postOfi = npost[i]
	for k in sparsevec[idxvec[i] : idxvec[i] + postOfi]:
	    # for angle in range(nPhis):
	    if k < NE:
		ue[k, :] += JEE * tc[i, :]
	    else:
		ue[k, :] += JIE * tc[i, :]
    plt.figure()
    poe = POofPopulation(ue, IF_IN_RANGE=1)
    plt.hist(poe[:NE], normed = 1, histtype = 'step')
    plt.hist(poe[NE:], normed = 1, histtype = 'step')
    return poe

def PlotPrePODistr(idxvec, sparsevec, npost, tc, JEE = 1, JEI = -1.5, JIE = 1, JII = -1, NE = 10000, NI = 10000, nPhis = 8):
    fr = tc.mean(1) # mean over angles
    po = POofPopulation(tc, IF_IN_RANGE = 1)
    # poTest = 180 * np.random.rand(NE)
    # po = poTest
    nNeurons = idxvec.size
    ue = np.zeros((nNeurons, nPhis))
    poOfPre = [[] for i in range(nNeurons)]
    for i in range(NE):
	postOfi = npost[i]
	for k in sparsevec[idxvec[i] : idxvec[i] + postOfi]:
	    poOfPre[k].append(po[i])
    plt.figure()
    plt.hist(np.hstack(poOfPre), normed = 1, histtype = 'step')
    fre = fr[:NE]
    frBiased = fre[po[:NE] > 120]
    plt.figure()
    plt.hist(fre, 23, normed = 1, histtype = 'step', label = 'full')
    plt.hist(frBiased, 23, normed = 1, histtype = 'step', label = 'biased')
    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 12})		
    return poOfPre, fr, frBiased

def PlotPODistrOfInput(N, K, trNo = 800, nPhis=8, kappa = 0, p = 0, gamma = 0, nPop = 2, rewireType = 'rand'):
    if(nPop == 2):
        ui = GetInputTuningCurves(p=p, gamma=gamma, rewireType = rewireType, kappa=0, inputType='EI', N=N, K=K, trNo=trNo, nPhis=nPhis)
        ue = GetInputTuningCurves(p=p, gamma=gamma, rewireType = rewireType, kappa=0, inputType='EE', N=N, K=K, trNo=trNo, nPhis=nPhis)
        poe = POofPopulation(ue, IF_IN_RANGE=1) * 180 / np.pi
        poi = POofPopulation(-ui, IF_IN_RANGE=1) * 180 / np.pi
        plt.figure()
        plt.hist(poi[:N], histtype='step', normed=1, color='r')            
        plt.hist(poe[:N], histtype='step', normed=1, color='k')
        legendtxt = [r'$u^i_{EE}(\phi)$', r'-1 x $u^i_{EI}(\phi)$']
        plt.title('N=%s, K = %s'%(N, K));
        plt.xlabel('PO (deg)')
        plt.legend(legendtxt, loc =8,  frameon = False, numpoints = 1, prop = {'size': 20}, markerscale=1, ncol=1)
        filename = 'po_uE_N%s_K%s'%(N, K)
        ProcessFigure(plt.gcf(), filename, 1, paperSize = [1.5*4, 1.5*3], axPosition=[.21, .15, .7, .65], titleSize=14, nDecimalsY=3, figFormat='eps')
        plt.figure()
        plt.hist(poe[N:], histtype='step', normed=1, color='k')
        plt.hist(poi[N:], histtype='step', normed=1, color='r')
        legendtxt = [r'$u^i_{IE}(\phi)$', r'-1 x $u^i_{II}(\phi)$']
        plt.title('N=%s, K = %s'%(N, K));
        plt.xlabel('PO (deg)') 
        legendtxt = [r'$u^i_{IE}(\phi)$', r'-1 x $u^i_{II}(\phi)$']   
        plt.legend(legendtxt, loc =8,  frameon = False, numpoints = 1, prop = {'size': 20}, markerscale=1, ncol=1)    
        filename = 'po_uI_N%s_K%s'%(N, K)
        ProcessFigure(plt.gcf(), filename, 1, paperSize = [1.5*4, 1.5*3], axPosition=[.21, .15, .7, .65], titleSize=14, nDecimalsY=3, figFormat='eps')
    else:
        ui = GetInputTuningCurves(p=p, gamma=gamma, kappa=0, inputType='EI', N=N, K=K, trNo=trNo, nPhis=nPhis)
        poi = POofPopulation(-ui, IF_IN_RANGE=1) * 180 / np.pi
        # plt.figure()
        plt.hist(poi[:N], histtype='step', normed=1, label = 'N=%s, K = %s'%(N, K))            
        plt.xlabel('PO (deg)') 
        titleTxt =  r'-1 x $u^i_{II}(\phi)$'
        plt.title(titleTxt);
        plt.legend(loc =8,  frameon = False, numpoints = 1, prop = {'size': 10}, markerscale=1, ncol=1)    
        filename = 'po_onepop_uI_N%s_K%s'%(N, K)
        ProcessFigure(plt.gcf(), filename, 1, paperSize = [1.5*4, 1.5*3], axPosition=[.21, .15, .7, .65], titleSize=14, nDecimalsY=3, figFormat='eps')

def CosFunc(x, a, b, c):
    #    x and c are in degrees
    x = x * np.pi / 180.0
    c = c * np.pi / 180.0
    return a + b * np.cos(2 * (x - c))

def PltFitHistCos(counts, bins, color):
    x = (bins[:-1] + bins[1:]) / 2.0
    nBins = x.size
    m0 = counts.mean()
    m1 = M1Component(counts) 
    phase = GetPhase(counts, x, IF_IN_RANGE = True)
    xx = np.arange(0, 180, 180./100)
    plt.plot(xx, CosFunc(xx, m0, m1, phase), '--', color = color)
    return m0, m1, phase

def PlotPODistrOfInputConducatance(ue, ui, N = 10000, K = 1000, trNo = 800, nPhis=8, kappa = 0, p = 0, gamma = 0, nPop = 2):
    if(nPop == 2):
        # ue = cc.GetTuningCurves('conductanceE_xi1.2_theta')
        # ui = cc.GetTuningCurves('conductanceE_xi1.2_theta')
        poe = POofPopulation(ue, IF_IN_RANGE=1)
        poi = POofPopulation(ui, IF_IN_RANGE=1)
        plt.figure()
        nCountEI, binsEI, _  = plt.hist(poi[:N], histtype='step', normed=1, color='r')            
        PltFitHistCos(nCountEI, binsEI, 'r')
        nCountEE, binsEE, _  =plt.hist(poe[:N], histtype='step', normed=1, color='k')
        PltFitHistCos(nCountEE, binsEE, 'k')
        legendtxt = [r'$g^i_{EE}(\phi)$', r'$g^i_{EI}(\phi)$']
        plt.title('N=%s, K = %s'%(N, K));
        plt.xlabel('PO (deg)')
        plt.legend(legendtxt, loc =8,  frameon = False, numpoints = 1, prop = {'size': 20}, markerscale=1, ncol=1)
        filename = 'po_uE_N%s_K%s_conductance'%(N, K)
        plt.ylim(0, plt.ylim()[1] + 0.002)
        ProcessFigure(plt.gcf(), filename, 1, paperSize = [1.5*4, 1.5*3], axPosition=[.21, .15, .7, .65], titleSize=14, nDecimalsX=3, figFormat='eps')
        plt.figure()
        nCountIE, binsIE, _  =plt.hist(poe[N:], histtype='step', normed=1, color='k')
        PltFitHistCos(nCountIE, binsIE, 'k')
        nCountII, binsII, _  =plt.hist(poi[N:], histtype='step', normed=1, color='r')
        PltFitHistCos(nCountII, binsII, 'r')

        legendtxt = [r'$g^i_{IE}(\phi)$', r'$g^i_{II}(\phi)$']
        plt.title('N=%s, K = %s'%(N, K));
        plt.xlabel('PO (deg)') 
        legendtxt = [r'$g^i_{IE}(\phi)$', r'$g^i_{II}(\phi)$']   
        plt.legend(legendtxt, loc =8,  frameon = False, numpoints = 1, prop = {'size': 20}, markerscale=1, ncol=1)    
        filename = 'po_uI_N%s_K%s_conductance'%(N, K)
        plt.ylim(0, plt.ylim()[1] + 0.002)
        ProcessFigure(plt.gcf(), filename, 1, paperSize = [1.5*4, 1.5*3], axPosition=[.21, .15, .7, .65], titleSize=14, nDecimalsX=3, figFormat='eps')

def PlotInputOSIDistr(ue, ui, N, K, filetag = ''):
    osiue = OSIOfPop(ue, 0)
    if np.any(ui < 0):
	ui = -ui
    osiui = OSIOfPop(ui, 0)
    osiE = osiue[:N]; osiE = osiE[~np.isnan(osiE)]
    osiI = osiui[:N]; osiI = osiI[~np.isnan(osiI)]
    print 'E', np.nanmean(osiE), np.nanmean(osiI)    
    plt.hist(osiE, 26,  normed = 1, histtype = 'step', label= 'EE', color = 'k')
    plt.hist(osiI, 25,  normed = 1, histtype = 'step', label= 'EI', color = 'r')
    plt.vlines(np.nanmean(osiue[:N]), *plt.ylim(), color = 'k')
    plt.vlines(np.nanmean(osiui[:N]), *plt.ylim(), color = 'r')
    filename = 'osi_uE_N%s_K%s_conductance'%(N, K) + filetag
    plt.xlabel('OSI')
    plt.xlim(0, .03)
    if np.nanmax(osiE) > 0.03:
	plt.xlim(0, 1)
    plt.legend(loc =0,  frameon = False, numpoints = 1, prop = {'size': 10}, markerscale=1, ncol=1)    
    plt.show()    
    ProcessFigure(plt.gcf(), filename, 1, paperSize = [4, 3], axPosition=[.21, .15, .7, .65], titleSize=14, nDecimalsX=3, figFormat='eps')
    # uI
    plt.clf()
    osiE = osiue[N:]; osiE = osiE[~np.isnan(osiE)]
    osiI = osiui[N:]; osiI = osiI[~np.isnan(osiI)]    
    plt.hist(osiE, 26,  normed = 1, histtype = 'step', label= 'IE', color = 'k')
    print 'I', np.nanmean(osiE), np.nanmean(osiui[N:])
    plt.hist(osiI,  25, normed = 1, histtype = 'step', label= 'II', color = 'r')
    plt.vlines(np.nanmean(osiue[N:]), *plt.ylim(), color = 'k')
    plt.vlines(np.nanmean(osiui[N:]), *plt.ylim(), color = 'r')
    filename = 'osi_uI_N%s_K%s_conductance'%(N, K) + filetag 
    plt.xlabel('OSI')
    plt.xlim(0, .03)
    if np.nanmax(osiE) > 0.03:
	plt.xlim(0, 1)
    plt.legend(loc =0,  frameon = False, numpoints = 1, prop = {'size': 10}, markerscale=1, ncol=1)    
    plt.show()    
    ProcessFigure(plt.gcf(), filename, 1, paperSize = [4, 3], axPosition=[.21, .15, .7, .65], titleSize=14, nDecimalsX=3, figFormat='eps')
    plt.clf()

def GetUEorIMinMat(x):
    #if x is ui then pass -ui
    N, nAngles = x.shape
    if np.any(x[:] < 0):
        x = -x
    xmin = np.min(x, 1)
    xmin.shape = N, 1
    z = np.tile(xmin, 8)
    xmin = x - z
    return xmin
    
def ComparePOofInOut( rewireType, p = 0, gamma = 0, mExt = 0.075, mExtOne = 0.075, nPhis = 8, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 0, IF_LEGEND = False):
    plt.figure()
    ui = GetInputTuningCurves(p=p, gamma=gamma, rewireType = rewireType, kappa=0, inputType='EI', N=N, K=K, trNo=trNo, nPhis=nPhis, T = T)
    ue = GetInputTuningCurves(p=p, gamma=gamma, rewireType = rewireType, kappa=0, inputType='EE', N=N, K=K, trNo=trNo, nPhis=nPhis, T = T)
    poe = POofPopulation(ue, IF_IN_RANGE = True)
    poi = POofPopulation(-ui, IF_IN_RANGE = True)
    tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
    po = POofPopulation(tc, IF_IN_RANGE = True)

    #
    pdf, bins, _ = plt.hist(poe[:N], normed = 1)
    pdfi, binsi, _ = plt.hist(poi[:N], normed = 1)
    plt.clf()
    binsE = (bins[:-1] + bins[1:]) / 2.0
    binsI = (binsi[:-1] + binsi[1:]) / 2.0

    #

    theta = np.arange(0, 180.0, 180.0 / nPhis)
    mEOftheta = tc[:N].mean(0)
    mIOftheta = tc[N:].mean(0)
    poOfME = GetPhase(mEOftheta, theta, IF_IN_RANGE=True)
    poOfMI = GetPhase(mIOftheta, theta, IF_IN_RANGE=True)
    poOfUEE = GetPhase(pdf, np.linspace(0, 180.0, pdf.size, endpoint = False), IF_IN_RANGE=True)
    poOfUEI = GetPhase(pdfi, np.linspace(0, 180.0, pdf.size, endpoint = False), IF_IN_RANGE=True)    
    plt.plot(binsE, pdf, 'k-', label  = r'$u_{EE}$')
    plt.plot(binsI, pdfi, 'r-', label =  r'$u_{EI}$')
    plt.vlines(poOfME, *plt.ylim(), lw = 2, color = 'k', linestyles='--', label = r'$PO(m_E)$')
    plt.vlines(poOfMI, *plt.ylim(), lw = 2, color = 'r', linestyles='--', label = r'$PO(m_I)$')    
    plt.vlines(poOfUEE, *plt.ylim(), lw = 2, color = 'k', label = r'$PO(u_{EE})$')
    plt.vlines(poOfUEI, *plt.ylim(), lw = 2, color = 'r', label = r'$PO(u_{EI})$')
    plt.title('realization #%s'%(trNo - 800 + 1))
    if IF_LEGEND:
        plt.legend(loc =0,  frameon = False, numpoints = 1, prop = {'size': 10}, markerscale=1, ncol=1)
    filename = 'ComparePOofInOut_N%s_K%s_tr%s_T%s'%(N, K, trNo, int(1e-3*T))
    ProcessFigure(plt.gcf(), filename, 1, paperSize = [4, 3], axPosition=[.21, .15, .7, .65], titleSize=14, nDecimalsX=0, nDecimalsY=3, figFormat='eps')    
    

def ComparePOofInOutScatter( rewireType, p = 0, gamma = 0, mExt = 0.075, mExtOne = 0.075, nPhis = 8, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 0, IF_LEGEND = False):
    plt.figure()
    poUE = []
    poUI = []
    poME = []
    poMI = []
    for trNo in range(800, 811):
        ui = GetInputTuningCurves(p=p, gamma=gamma, rewireType = rewireType, kappa=0, inputType='EI', N=N, K=K, trNo=trNo, nPhis=nPhis, T = T)
        ue = GetInputTuningCurves(p=p, gamma=gamma, rewireType = rewireType, kappa=0, inputType='EE', N=N, K=K, trNo=trNo, nPhis=nPhis, T = T)
        poe = POofPopulation(ue, IF_IN_RANGE = True)
        poi = POofPopulation(-ui, IF_IN_RANGE = True)
        tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
        po = POofPopulation(tc, IF_IN_RANGE = True)

        #
        pdf, bins, _ = plt.hist(poe[:N], normed = 1)
        pdfi, binsi, _ = plt.hist(poi[:N], normed = 1)
        plt.clf()
        binsE = (bins[:-1] + bins[1:]) / 2.0
        binsI = (binsi[:-1] + binsi[1:]) / 2.0

        #

        theta = np.arange(0, 180.0, 180.0 / nPhis)
        mEOftheta = tc[:N].mean(0)
        mIOftheta = tc[N:].mean(0)
        poME.append(GetPhase(mEOftheta, theta, IF_IN_RANGE=True))
        poMI.append(GetPhase(mIOftheta, theta, IF_IN_RANGE=True))
        poUE.append(GetPhase(pdf, np.linspace(0, 180.0, pdf.size, endpoint = False), IF_IN_RANGE=True))
        poUI.append(GetPhase(pdfi, np.linspace(0, 180.0, pdf.size, endpoint = False), IF_IN_RANGE=True))
        
    # if IF_LEGEND:
    #     plt.legend(loc =0,  frameon = False, numpoints = 1, prop = {'size': 10}, markerscale=1, ncol=1)

    x = np.arange(0, 180, 180.0 / 100)
    plt.plot(poME, poUE, 'ko')
    plt.xlabel('Peak of ' + r'$m_E(\theta)$')
    plt.ylabel('Peak of ' + r'$PO(u_{EE})$' + ' pdf')    
    filename = 'Scatter_E_ComparePOofInOut_N%s_K%s_tr%s_T%s'%(N, K, trNo, int(1e-3*T))
    plt.plot(x, x, 'k-')
    ProcessFigure(plt.gcf(), filename, 1, paperSize = [4, 3], axPosition=[.22, .2, .66, .65], titleSize=14, nDecimalsX=0, nDecimalsY=0, figFormat='eps')
    plt.figure()
    plt.plot(poMI, poUI, 'ro')
    plt.xlabel('Peak of ' + r'$m_I(\theta)$')
    plt.ylabel('Peak of ' + r'$PO(u_{EI})$' + ' pdf')
    plt.plot(x, x, 'k-')
    filename = 'Scatter_I_ComparePOofInOut_N%s_K%s_tr%s_T%s'%(N, K, trNo, int(1e-3*T))
    ProcessFigure(plt.gcf(), filename, 1, paperSize = [4, 3], axPosition=[.22, .2, .65, .65], titleSize=14, nDecimalsX=0, nDecimalsY=0, figFormat='eps')

    
def M1vsSimTimeChks(nChunks, rewireType, p = 0, gamma = 0, nPhis = 8, mExt = 0.075, mExtOne = 0.075, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 0, IF_SUCCESS = False):
    m1EofAvg = np.empty((nChunks, ))
    m1IofAvg = np.empty((nChunks, ))
    m1EofAvg[:] = np.nan; m1IofAvg[:] = np.nan;
    M1Func = lambda x, theta: 2.0 * np.absolute(np.dot(x, np.exp(-2.0j * theta))) / len(x)
    theta = np.arange(0, np.pi, np.pi / nPhis)
    simT = np.linspace(0, 1, nChunks + 1, endpoint = False)
    mEOld = np.zeros((nPhis, )); mIOld = np.zeros((nPhis, ))
    m0EofrunAvg = np.zeros((nPhis, ))
    m0IofrunAvg = np.zeros((nPhis, ))
    for i in range(nChunks):
        tc = GetTuningCurvesChnk(i, rewireType, p, gamma, nPhis, mExt, mExtOne, trNo, N, K, nPop, T, kappa, IF_SUCCESS)
        if i == 0:
            m0EofrunAvg = np.nanmean(tc[:N], 0)
            m0IofrunAvg = np.nanmean(tc[N:], 0)
        else:
            m0EofrunAvg = 0.5 * (m0EofrunAvg + np.nanmean(tc[:N], 0))
            m0IofrunAvg = 0.5 * (m0IofrunAvg + np.nanmean(tc[N:], 0))            
        vidx = ~np.isnan(m0EofrunAvg)
        m1EofAvg[i] = M1Func(m0EofrunAvg[vidx], theta[vidx])
        m1IofAvg[i] = M1Func(m0IofrunAvg[vidx], theta[vidx])
    xaxis = 1.0 / np.sqrt(simT[1:])
    plt.plot(simT[:-1], m1EofAvg, 'ks-', label = 'E')
    plt.plot(simT[:-1], m1IofAvg, 'rs-', label = 'I')
    plt.xlabel(r'$\frac{1}{\sqrt{T}}$')
    plt.ylabel('modulation')

    
def M1vsSimTime(rewireType, p = 0, gamma = 0, nPhis = 8, mExt = 0.075, mExtOne = 0.075, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 0, IF_SUCCESS = False):
    TList = [1000, 2000, 4000, 6000, 8000, 10000]
    nSimTimes = len(TList)
    nChunks = nSimTimes
    m1EofAvg = np.empty((nChunks, ))
    m1IofAvg = np.empty((nChunks, ))
    m1EofAvg[:] = np.nan; m1IofAvg[:] = np.nan;
    M1Func = lambda x, theta: 2.0 * np.absolute(np.dot(x, np.exp(-2.0j * theta))) / len(x)
    theta = np.arange(0, np.pi, np.pi / nPhis)
    simT = np.linspace(0, 1, nChunks + 1, endpoint = False)
    mEOld = np.zeros((nPhis, )); mIOld = np.zeros((nPhis, ))
    m0EofrunAvg = np.zeros((nPhis, ))
    m0IofrunAvg = np.zeros((nPhis, ))
    fitFuncLint = lambda x, m, c: m * x  + c
    for i, T in enumerate(TList):
        # tc = GetTuningCurvesChnk(i, rewireType, p, gamma, nPhis, mExt, mExtOne, trNo, N, K, nPop, T, kappa, IF_SUCCESS)
        tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)        
        if i == 0:
            m0EofrunAvg = np.nanmean(tc[:N], 0)
            m0IofrunAvg = np.nanmean(tc[N:], 0)
        else:
            m0EofrunAvg = 0.5 * (m0EofrunAvg + np.nanmean(tc[:N], 0))
            m0IofrunAvg = 0.5 * (m0IofrunAvg + np.nanmean(tc[N:], 0))            
        vidx = ~np.isnan(m0EofrunAvg)
        m1EofAvg[i] = M1Func(m0EofrunAvg[vidx], theta[vidx])
        m1IofAvg[i] = M1Func(m0IofrunAvg[vidx], theta[vidx])

    x = 1.0 / np.sqrt(np.array(TList))
    plt.plot(x, m1EofAvg, 'ks-', label = r'$m_E^{(1)}$')
    plt.plot(x, m1IofAvg, 'rs-', label = r'$m_I^{(1)}$')
    slope0 = 1; intercept0 = 0.1
    fitParams, fitError = curve_fit(fitFuncLint, x, m1EofAvg, p0 = [slope0, intercept0])
    xx = np.linspace(0, x[0], 100)
    print 'E:', fitParams
    plt.plot(xx, fitFuncLint(xx, *fitParams), 'k--')
    fitParams, fitError = curve_fit(fitFuncLint, x, m1IofAvg, p0 = [slope0, intercept0]) 
    print 'I:', fitParams    
    plt.plot(xx, fitFuncLint(xx, *fitParams), 'r--')


    # plt.plot(np.sqrt(np.array(TList)), m1EofAvg, 'ks-', label = 'E')
    # plt.plot(np.sqrt(np.array(TList)), m1IofAvg, 'rs-', label = 'I')

    plt.xlabel(r'$\frac{1}{\sqrt{T}}$')
    plt.ylabel(r'$\mathrm{modulation\, of\,  } \left\langle m_i  \right\rangle_i (\theta)$', fontsize = 16)
    # plt.ylim(0, plt.ylim()[1])
    plt.legend(loc =0,  frameon = False, numpoints = 1, prop = {'size': 10}, markerscale=1, ncol=1)        
    filename = 'avgActM1_vs_T_N%s_K%s_tr%s'%(N, K, trNo)
    ProcessFigure(plt.gcf(), filename, 1, paperSize = [4, 3], axPosition=[.29, .26, .65, .65], titleSize=14, nDecimalsX=3, nDecimalsY=4, figFormat='eps', labelFontsize = 16)


def M1ofUEvsSimTime(rewireType, p = 0, gamma = 0, nPhis = 8, mExt = 0.075, mExtOne = 0.075, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 0, IF_SUCCESS = False):
    TList = [1000, 2000, 4000, 6000, 8000, 10000]
    nSimTimes = len(TList)
    nChunks = nSimTimes
    m1EofAvg = np.empty((nChunks, ))
    m1IofAvg = np.empty((nChunks, ))
    m1EofAvg[:] = np.nan; m1IofAvg[:] = np.nan;
    M1Func = lambda x, theta: 2.0 * np.absolute(np.dot(x, np.exp(-2.0j * theta))) / len(x)
    theta = np.arange(0, np.pi, np.pi / nPhis)
    simT = np.linspace(0, 1, nChunks + 1, endpoint = False)
    mEOld = np.zeros((nPhis, )); mIOld = np.zeros((nPhis, ))
    m0EofrunAvg = np.zeros((nPhis, ))
    m0IofrunAvg = np.zeros((nPhis, ))
    fitFuncLint = lambda x, m, c: m * x  + c

    m1ofFF = []
    m1ofUE = []

    cosFunc = lambda theta, m0, m1, phase: m0 + m1 * np.cos(2*(theta + phase))
    for i, T in enumerate(TList):
        try:
            uff = GetInputTuningCurves(p=p, gamma=gamma, rewireType = rewireType, kappa=0, inputType='FF', N=N, K=K, trNo=trNo, nPhis=nPhis, T = T)
            ue = GetInputTuningCurves(p=p, gamma=gamma, rewireType = rewireType, kappa=0, inputType='EE', N=N, K=K, trNo=trNo, nPhis=nPhis, T = T)
            poff = POofPopulation(uff, IF_IN_RANGE = True)            
            poe = POofPopulation(ue, IF_IN_RANGE = True)
            try:
                pdf, bins, _ = plt.hist(poe[:N], normed = 1, histtype = 'step')
                pdfFF, binsFF, _ = plt.hist(poff[:N], normed = 1, histtype = 'step')
                
                # plt.hist(poff[:N], normed = 1)                
                fitParams, fitError = curve_fit(cosFunc, bins[:-1]*np.pi / 180, pdf, p0 = [.01, .001, np.pi / 2.0])                
                m1ofUE.append(fitParams[1])
                fitParams, fitError = curve_fit(cosFunc, binsFF[:-1]*np.pi / 180, pdfFF, p0 = [.01, .001, np.pi / 2.0])
                m1ofFF.append(fitParams[1])
                
                PltFitHistCos(pdfFF, binsFF, 'r')
                xx = np.linspace(0, np.pi, 100)
                plt.plot(xx * 180.0 / np.pi, cosFunc(xx, *fitParams), 'c--')
                print 'fit:', fitParams
                print 'est:', pdfFF.mean(), M1Component(pdfFF)
                # plt.title('%s'%(T))
                plt.ylim(0, .009)
                # plt.waitforbuttonpress()
                plt.clf()
            except AttributeError:
                m1ofUE.append(np.nan); m1ofFF.append(np.nan)
        except IOError:
            m1ofUE.append(np.nan); m1ofFF.append(np.nan)
    x = 1.0 / np.sqrt(np.array(TList))
    plt.plot(x, m1ofUE, 'ks-', label = 'EE')
    plt.plot(x, m1ofFF, 'gs-', label = 'FF')
    slope0 = 1; intercept0 = 0.1
    fitParams, fitError = curve_fit(fitFuncLint, x, m1ofUE, p0 = [slope0, intercept0])
    xx = np.linspace(0, x[0], 100)
    print ' E:', fitParams
    plt.plot(xx, fitFuncLint(xx, *fitParams), 'k--')
    fitParams, fitError = curve_fit(fitFuncLint, x, m1ofFF, p0 = [slope0, intercept0]) 
    print 'FF:', fitParams    
    plt.plot(xx, fitFuncLint(xx, *fitParams), 'g--')
    plt.xlabel(r'$\frac{1}{\sqrt{T}}$')
    # plt.ylabel(r'$\mathrm{modulation\, of\,  } \left\langle m_i  \right\rangle_i (\theta)$', fontsize = 16)
    plt.ylabel(r'modulation of PO distr', fontsize = 10)
    plt.legend(loc =0,  frameon = False, numpoints = 1, prop = {'size': 10}, markerscale=1, ncol=1)    
    # plt.ylim(0, plt.ylim()[1])
    
    filename = 'POdistr_M1_vs_T_N%s_K%s_tr%s'%(N, K, trNo)
    ProcessFigure(plt.gcf(), filename, 1, paperSize = [4, 3], axPosition=[.26, .26, .65, .65], titleSize=14, nDecimalsX=3, nDecimalsY=4, figFormat='eps', labelFontsize = 12)    

    

def M1vsKoverN(rewireType, p = 0, gamma = 0, nPhis = 8, mExt = 0.075, mExtOne = 0.075, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 0, IF_SUCCESS = False):
    # KList = np.array([500, 1000, 2000, 3000, 4000])
    # KList = np.array([500])
    NList = np.array([5000, 10000, 20000, 30000, 40000])
    # NList = np.array([40000])    
    # KList = np.array([2000])
    # KoverNList = KList / float(N)
    m1ofFF = []
    m1ofUE = []
    # for i, K in enumerate(KList):
    for i, N in enumerate(NList):        
        try:
            tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
            uff = GetInputTuningCurves(p=p, gamma=gamma, rewireType = rewireType, kappa=0, inputType='FF', N=N, K=K, trNo=trNo, nPhis=nPhis, T = T)

            ue = GetInputTuningCurves(p=p, gamma=gamma, rewireType = rewireType, kappa=0, inputType='EE', N=N, K=K, trNo=trNo, nPhis=nPhis, T = T)
            poff = POofPopulation(uff, IF_IN_RANGE = True)            
            poe = POofPopulation(ue, IF_IN_RANGE = True)
            try:
                pdf, bins, _ = plt.hist(poe[:N], normed = 1)
                pdfFF, binsFF, _ = plt.hist(poff[:N], normed = 1)
                m1ofFF.append(M1Component(pdfFF))
                m1ofUE.append(M1Component(pdf))
                # plt.waitforbuttonpress()
                plt.clf()
            except AttributeError:
                m1ofUE.append(np.nan); m1ofFF.append(np.nan)
        except IOError:
            m1ofUE.append(np.nan); m1ofFF.append(np.nan)
    print m1ofUE
    x = 1.0 / np.sqrt(NList)
    plt.plot(x, m1ofFF, 'ko--', label = 'FF')
    plt.plot(x, m1ofUE, 'ks-', label = 'EE')
    # plt.plot(KoverNList, m1ofFF, 'ko--', label = 'FF')
    # plt.plot(KoverNList, m1ofUE, 'ks-', label = 'EE')
    plt.legend(loc =0,  frameon = False, numpoints = 1, prop = {'size': 10}, markerscale=1, ncol=1)
    plt.xlabel(r'$\frac{1}{\sqrt{N}}$')    
    plt.ylabel('modulation of PO distr')
    filename = 'POdistr_M1_vs_N_T%s_K%s_tr%s'%(int(T*1e-3), K, trNo)
    ProcessFigure(plt.gcf(), filename, 1, paperSize = [4, 3], axPosition=[.26, .26, .65, .65], titleSize=14, nDecimalsX=3, nDecimalsY=4, figFormat='eps', labelFontsize = 12)    


def RandomSum(tce, filetag, IF_PLOT, K):
    nNeurons, nPhis = tce.shape
    osi = OSIOfPop(tce, 0)
    m0 = tce.mean();
    m1 = 2.0 * np.nanmean(osi) * m0
    theta = np.arange(0, np.pi, np.pi / 8.0)
    tcControl = lambda neuronIdx : m0 + m1 * np.cos(2.0 * (theta - neuronIdx * np.pi / nNeurons))
    rndInputsCntrl = np.zeros((nNeurons, nPhis))    
    rndInputs = np.zeros((nNeurons, nPhis))
    neuronIdx = range(nNeurons)
    m1ofControl = []; m1ofRandSum = []
    m1ofPOdistrControl = []; m1ofPOdistrRandSum = []
    # ipdb.set_trace()
    for i in range(nNeurons):
        rndPreNeurons = np.random.choice(neuronIdx, K, replace = False)
        for k in rndPreNeurons:
            rndInputsCntrl[i, :] += tcControl(k)
        rndInputs[i, :] += np.sum(tce[rndPreNeurons, :], 0)
    poCntrl = POofPopulation(rndInputsCntrl, IF_IN_RANGE = 1)
    po = POofPopulation(rndInputs, IF_IN_RANGE = 1)
    # ipdb.set_trace()    
    try:
        plt.ioff()
        pdfFF, binsFF, _ = plt.hist(po, normed = 1)
        pdf, bins, _ = plt.hist(poCntrl, normed = 1)   
        fitParams, fitError = curve_fit(CosFunc, bins[:-1]*np.pi / 180, pdf, p0 = [.01, .001, np.pi / 2.0])                
        m1ofPOdistrControl = M1Component(pdf) #fitParams[1]
        fitParams, fitError = curve_fit(CosFunc, binsFF[:-1]*np.pi / 180, pdfFF, p0 = [.01, .001, np.pi / 2.0])
        m1ofPOdistrRandSum = M1Component(pdfFF)#fitParams[1]
        m1ofRandSum = np.nanmean(PopM1Component(rndInputs)) #(pdf)
        m1ofControl = np.nanmean(PopM1Component(rndInputsCntrl))  #M1Component(pdfFF)
	# phase = GetPhase(pdfFF, bins[:-1])
	# theta = np.linspace(0, 180, 100)
	# plt.plot(theta, CosFunc(theta, pdfFF.mean(), M1Component(pdfFF), phase))
    except AttributeError:
        m1ofRandSum.append(np.nan); m1ofControl.append(np.nan)
    # ipdb.set_trace()
    IF_PLOT = False
    if IF_PLOT:
        plt.ioff()
        plt.clf()
        plt.hist(poCntrl, normed = 1, histtype = 'step', label = 'control')    
        plt.hist(po, normed = 1, histtype = 'step', label = 'simulation')
        plt.legend(loc =0,  frameon = False, numpoints = 1, prop = {'size': 10}, markerscale=1, ncol=1)

        plt.xlabel(r'$\mathrm{PO}$')    
        plt.ylabel('PDF')
        filename = 'rndsum_uE_N%s_K%s_tr%s'%(nNeurons, K, filetag)
        plt.title('realization#%s, K = %s'%(filetag, K))
        ProcessFigure(plt.gcf(), filename, 1, paperSize = [4, 3], axPosition=[.22, .22, .65, .65], titleSize=14, nDecimalsX=0, nDecimalsY=3, figFormat='eps', labelFontsize = 12)
    return [m1ofControl, m1ofRandSum, m1ofPOdistrControl, m1ofPOdistrRandSum]

def RandomSumvsK(tce, KList = np.arange(500, 5000, 500), filetag = '', IF_COMPUTE = False):
    nNeurons, nPhis = tce.shape
    if IF_COMPUTE:
	m1ofControl = np.zeros((len(KList), )); m1ofControl[:] = np.nan
	m1ofRandSum = np.zeros((len(KList), )); m1ofRandSum[:] = np.nan
	m1ofPOdistrRandSum = np.zeros((len(KList), )); m1ofPOdistrRandSum[:] = np.nan
	m1ofPOdistrControl = np.zeros((len(KList), )); m1ofPOdistrControl[:] = np.nan	
	pool = Pool(len(KList))
	func = partial(RandomSum, tce, '', True)
	out = pool.map(func, KList)
	pool.close()
	# ipdb.set_trace()
	for i, K in enumerate(KList):
	    m1ofControl[i] = out[i][0]
	    m1ofRandSum[i] = out[i][1]
	    m1ofPOdistrControl[i] = out[i][2]
            m1ofPOdistrRandSum[i] = out[i][3]	    
	x = 1.0 / np.sqrt(np.array(KList))
	np.save('./data/twopop/pobias/randomSum_vs_K_N%s_trNo'%(nNeurons) + filetag, [x, m1ofControl, m1ofRandSum, m1ofPOdistrControl, m1ofPOdistrRandSum])
    else:
	x, m1ofControl, m1ofRandSum, m1ofPOdistrControl, m1ofPOdistrRandSum = np.load('./data/twopop/pobias/randomSum_vs_K_N%s_trNo%s.npy'%(nNeurons, filetag))
	xx = np.linspace(1 / np.sqrt(float(KList[-1])), 1 / np.sqrt(float(KList[0])), 10)

	slope0 = 1; intercept0 = 0.1
	fitFuncLint = lambda x, m, c: m * x  + c    
	fitParamsControl, fitError = curve_fit(fitFuncLint, x, m1ofControl, p0 = [slope0, intercept0])
	fitParamsRandSum, fitError = curve_fit(fitFuncLint, x, m1ofRandSum, p0 = [slope0, intercept0])

	
	plt.plot(KList, m1ofControl, 'k', label = 'control')
	plt.plot(KList, m1ofRandSum, 'g', label = 'rand sum')
	plt.legend(loc =0,  frameon = False, numpoints = 1, prop = {'size': 10}, markerscale=1, ncol=1)        
        plt.xlabel(r'$K$')                
        plt.ylabel(r'$\left\langle u_{EE, i}^{(1)} \right\rangle_i$')
	filename = './figs/twopop/pobias/rndsum_cntrl_K_vs_m1_N%s_tr%s'%(nNeurons, filetag)
	ProcessFigure(plt.gcf(), filename, 1, paperSize = [4, 3], axPosition=[.22, .22, .65, .65], titleSize=14, nDecimalsX=3, nDecimalsY=5, figFormat='eps', labelFontsize = 12)

        
        plt.figure()
	plt.plot(KList, m1ofPOdistrControl, 'k', label = 'control')
	plt.plot(KList, m1ofPOdistrRandSum, 'g', label = 'rand sum')
        plt.xlabel(r'$K$')
	plt.legend(loc =0,  frameon = False, numpoints = 1, prop = {'size': 10}, markerscale=1, ncol=1)
	plt.title('modulation of PO distr')	        
	filename = './figs/twopop/pobias/rndsum_cntrl_Kvs_POdistrm1_N%s_tr%s'%(nNeurons, filetag)
	ProcessFigure(plt.gcf(), filename, 1, paperSize = [4, 3], axPosition=[.22, .22, .65, .65], titleSize=14, nDecimalsX=0, nDecimalsY=4, figFormat='eps', labelFontsize = 12)
        
        plt.figure()
        y0 = m1ofPOdistrControl /  np.asarray(KList, dtype = float)
        y1 = m1ofPOdistrRandSum /  np.asarray(KList, dtype = float)
	plt.plot(x, y0, 'k', label = 'control')
	plt.plot(x, y1, 'g', label = 'rand sum')
	fitParamsControl, fitError = curve_fit(fitFuncLint, x, y0, p0 = [slope0, intercept0])
	fitParamsRandSum, fitError = curve_fit(fitFuncLint, x, y1, p0 = [slope0, intercept0])  
	plt.plot(xx, fitFuncLint(xx, *fitParamsControl), 'k--')    
	plt.plot(xx, fitFuncLint(xx, *fitParamsRandSum), 'g--')
        plt.xlabel(r'$\frac{1}{\sqrt{K}}$')        
	plt.title('modulation of PO distr(normalized)')	
	plt.legend(loc =0,  frameon = False, numpoints = 1, prop = {'size': 10}, markerscale=1, ncol=1)    
	print ' control:', fitParamsControl
	print 'rand sum:', fitParamsRandSum
	filename = './figs/twopop/pobias/rndsum_cntrl_N%s_tr%s'%(nNeurons, filetag)
	ProcessFigure(plt.gcf(), filename, 1, paperSize = [4, 3], axPosition=[.22, .22, .65, .65], titleSize=14, nDecimalsX=3, nDecimalsY=6, figFormat='eps', labelFontsize = 12)
    
def RMSOSI(tce, filetag, K):
    nNeurons, nPhis = tce.shape
    osi = OSIOfPop(tce, 0)
    m0 = tce.mean();
    m1 = 2.0 * np.nanmean(osi) * m0
    theta = np.arange(0, np.pi, np.pi / 8.0)
    tcControl = lambda neuronIdx : m0 + m1 * np.cos(2.0 * (theta - neuronIdx * np.pi / nNeurons))
    rndInputsCntrl = np.zeros((nNeurons, nPhis))    
    rndInputs = np.zeros((nNeurons, nPhis))
    neuronIdx = range(nNeurons)
    for i in range(nNeurons):
        rndPreNeurons = np.random.choice(neuronIdx, K, replace = False)
        rndInputs[i, :] += np.sum(tce[rndPreNeurons, :], 0)
    po = POofPopulation(rndInputs, IF_IN_RANGE = 1)
    uEE1 = PopM1Component(rndInputs)
    rmsUEE1 = np.sqrt(np.mean(uEE1**2)) / np.mean(uEE1)
    return [uEE1, rmsUEE1, po]

def RMSOSIvsK(tce, filetag, KList, IF_COMPUTE = False):
    nNeurons, nPhis = tce.shape
    if IF_COMPUTE:
        uEE1ofRandSum = np.zeros((len(KList), nNeurons));
        uEE1ofRandSum[:] = np.nan
        POuEE1ofRandSum = np.zeros((len(KList), nNeurons));
        POuEE1ofRandSum[:] = np.nan
        rmsUEE1ofRandSum = np.zeros((len(KList), ));
        rmsUEE1ofRandSum[:] = np.nan
        pool = Pool(len(KList))
        func = partial(RMSOSI, tce, '')
        out = pool.map(func, KList)
        pool.close()
        for i, K in enumerate(KList):
            uEE1ofRandSum[i, :] = out[i][0]
            rmsUEE1ofRandSum[i] = out[i][1]
            POuEE1ofRandSum[i, :] = out[i][0]            
        np.savez('./data/twopop/pobias/uEE1_RMS_vs_K_N%s_trNo'%(nNeurons) + filetag, klist = KList, uEE1ofRandSum = uEE1ofRandSum, rmsUEE1ofRandSum = rmsUEE1ofRandSum, POuEE1ofRandSum = POuEE1ofRandSum)
    else:
        out = np.load('./data/twopop/pobias/uEE1_RMS_vs_K_N%s_trNo'%(nNeurons) + filetag + '.npz')
        KList = out['klist']
        uEE1ofRandSum = out['uEE1ofRandSum']
        rmsUEE1ofRandSum = out['rmsUEE1ofRandSum']
        POuEE1ofRandSum = out['POuEE1ofRandSum']
        plt.plot(KList, rmsUEE1ofRandSum)
        plt.xlabel('K')
        plt.ylabel(r'$\mathrm{RMS} (u^{(1)}) / [u^{(1)}]$')

        nBins = 9
        poBins = np.linspace(0, 180, nBins + 1)
        binIndex = np.digitize(POuEE1ofRandSum[-1] * 180 / np.pi, poBins)
        meanOSIinBin = []
        uEE1ofLast = uEE1ofRandSum[-1]

        for i in range(nBins):
            meanOSIinBin.append(np.nanmean(uEE1ofLast[binIndex == i]))
        plt.figure()
        plt.plot(np.arange(0, 180, 180.0/nBins), meanOSIinBin)
        plt.xlabel('PO')
        plt.ylabel('mean OSI')
        
def PObiasBootStrapAux(tce, K, nBins, filename, poolIdx):
    nNeurons, nPhis = tce.shape
    rndInputs = np.zeros((nNeurons, nPhis))
    neuronIdx = range(nNeurons)
    np.random.seed(1234 + 12345 * poolIdx)
    # plusOrMinus = -1 if np.random.rand() > 0.5 else 1
    # Kpre = K + plusOrMinus * np.sqrt(K) * np.random.randn()
    Kpre = K
    print poolIdx
    for k in range(nNeurons):
        tce[k, :] = np.roll(tce[k, :], np.random.randint(0, nPhis))
    for i in range(nNeurons):
        rndPreNeurons = np.random.choice(neuronIdx, Kpre, replace = False)
        rndInputs[i, :] += np.sum(tce[rndPreNeurons, :], 0)
    poOfUEE = POofPopulation(rndInputs, IF_IN_RANGE = 1)
    pdfPOofUEE, bins = np.histogram(poOfUEE, nBins, density = True)
    m1OfPOdistr = M1Component(pdfPOofUEE)
    fp = open(filename, 'a')
    np.savetxt(fp, [m1OfPOdistr])
    fp.close()
    return m1OfPOdistr

def PObiasBootStrap(tce, ue, filetag, K = 1000, nBins = 10, nIterations = 100, IF_COMPUTE = 0):
    nNeurons, nPhis = tce.shape
    filename =  './data/twopop/pobias/bootstrap_pobias_N%s_nIters%s_trNo%s.csv'%(nNeurons, nIterations, filetag)
    if IF_COMPUTE:
        pool = Pool(48)
        func = partial(PObiasBootStrapAux, tce, K, nBins, filename)
        out = pool.map(func, np.arange(nIterations))
        pool.close()
        out = np.asarray(out, dtype = float)
        np.save('./data/twopop/pobias/bootstrap_pobias_N%s_nIters%s_trNo%s'%(nNeurons, nIterations, filetag), out)
    else:
        # out = np.load('./data/twopop/pobias/bootstrap_pobias_N%s_nIters%s_trNo%s.npy'%(nNeurons, nIterations, filetag))
        fp = open(filename)
        out = np.loadtxt(fp)
        fp.close()
        # ipdb.set_trace()
        plt.hist(out, 20, normed = 1, histtype = 'step', color = 'k')
        poOfUEESim = POofPopulation(ue, IF_IN_RANGE = 1)
        pdfPOofUEESim, bins = np.histogram(poOfUEESim, nBins, density = True)
        m1OfPOdistrSim = M1Component(pdfPOofUEESim)
        # plt.vlines(m1OfPOdistrSim, *plt.ylim(), lw = 2, color = 'r')
        # confidence interval
        stdErrorSamplingDistr = 2 * np.std(out)
        print stdErrorSamplingDistr
        # plt.text(0.006, 200, 'CI = %.6s'%(stdErrorSamplingDistr), fontsize = 16)
        print 'mean=', np.mean(out)
        plt.vlines(np.mean(out), *plt.ylim(), lw = 2, color = 'r', label = 'mean')
	plt.legend(loc =0,  frameon = False, numpoints = 1, prop = {'size': 16}, markerscale=1, ncol=1)
        plt.xlabel('F1 of PO distr')
        plt.ylabel('Density')
        plt.title('F1 sampling distribution')
        plt.savefig('F1_sampling_distr.png')
        plt.savefig('F1_sampling_distr.eps')        
        
def BootStrapSummary(bootstrapCI, nBins = 20):
    m1OfPOdistrSim = []
    for trNo in range(800, 811):
        ue = GetInputTuningCurves(0, rewireType='cntrl', kappa=0, inputType='EE', N=10000, K=1000, trNo=trNo, nPhis=8)
        poOfUEESim = POofPopulation(ue, IF_IN_RANGE = 1)
        pdfPOofUEESim, bins = np.histogram(poOfUEESim, nBins, density = True)
        m1OfPOdistrSim.append(M1Component(pdfPOofUEESim))
    plt.plot(np.zeros((len(m1OfPOdistrSim), )), m1OfPOdistrSim, 'o', markersize = 4, color = [0.65, .65, .65], label = 'simulation')
    (_, caps, _) = plt.errorbar(0, np.nanmean(m1OfPOdistrSim), yerr = bootstrapCI, color = 'k', fmt = 'o', markersize = 4, lw = 1, capsize = 20, elinewidth = 1, label = '   2 SD')
    for cap in caps:
        cap.set_markeredgewidth(.95)
    plt.xlim(-.06, .3)
    plt.legend(loc = 0,  frameon = False, numpoints = 1, prop = {'size': 10}, markerscale=.75, ncol=1)
    plt.xlabel('')
    plt.ylabel('F1 of PO distr')
    plt.show()
    plt.gca().set_xticklabels('')    
    filename = 'bootstrap_summary'
    ProcessFigure(plt.gcf(), filename, 1, paperSize = [4*.75, .75*3], axPosition=[.3, .22, .65, .65], titleSize=14, nDecimalsX=3, nDecimalsY=4, figFormat='eps', labelFontsize = 10)



def WritePOToTr0FolderBeforeRewiring(kappa, p, gamma, mExt, mExtOne, nPhis=8, rewireType='rand', N=10000, K=1000, nPop=2, T=1000, trNo = 0):
    out = np.nan
    try:
	tcOut = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
        po = POofPopulation(tcOut, IF_IN_RANGE = 1)
        requires = ['CONTIGUOUS', 'ALIGNED']
        po = np.require(po, np.float64, requires) * np.pi / 180.0
        baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
        fp = open(baseFldr + 'poOfNeurons.dat', 'wb')
        po.tofile(fp)
        fp.close()
        # for ii in range(100):
        #     theta = np.arange(0, 180, 22.5)
        #     plt.plot(theta, tcOut[ii, :])
        #     plt.vlines(po[ii] * 180 / np.pi, *plt.ylim())
        #     plt.waitforbuttonpress()
        #     plt.clf()
        print po[:11]
    except IOError:
        print 'file not found'
    

    
        
