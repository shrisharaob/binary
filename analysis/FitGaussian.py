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
from scipy.special import gammaincc as GammaQAux
from scipy.signal import argrelextrema
import Scripts as sc
from lmfit import minimize, Parameters


def GammaQ(a, x):
    return GammaQAux(0.5 * a, 0.5 * x)

def LoadFrChnk(p, gamma, phi, mExt, mExtOne, chkId, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False):
    baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
    if nPop == 1:
    	baseFldr = baseFldr + 'onepop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if nPop == 2:
    	baseFldr = baseFldr + 'twopop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if IF_VERBOSE:
    	print baseFldr
    filename = 'meanrates_theta%.6f_tr%s_chnk%s.txt'%(phi, trNo, chkId)
    return np.loadtxt(baseFldr + filename)

def LoadFr(p, gamma, phi, mExt, mExtOne, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False):
    baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
    if nPop == 1:
    	baseFldr = baseFldr + 'onepop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if nPop == 2:
    	baseFldr = baseFldr + 'twopop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if IF_VERBOSE:
    	print baseFldr
    filename = 'meanrates_theta%.6f_tr%s_last.txt'%(phi, trNo)
    return np.loadtxt(baseFldr + filename)

def Residual(params, theta, data, dataSEM):
    amp = params['amp']
    po = params['po']
    sigma = params['sigma']
    offset = params['offset']
    po = AngleInRange(po)
    sigma = AngleInRange(sigma)
    model = 0.0
    for n in np.arange(-2, 2):
	model += amp * np.exp(-(theta - po + 180.0 * n)**2 / (sigma**2)) + offset
    return (data - model) #/ dataSEM

def GetParams(params):
    amp = params['amp'].value
    po = params['po'].value
    sigma = params['sigma'].value
    offset = params['offset'].value
    return AngleInRange(po), AngleInRange(sigma), amp, offset

def FitPeriodicGaussianNEW(x, y, ySEM, p0 = [90, 90, 0.1, 1]):
    params = Parameters()
    # params.add('po', value = p0[0], min = 0.0, max = 180.0)
    # params.add('sigma', value = p0[1], min = 0.0, max = 180.0)
    # params.add('amp', value = p0[2], min = 0.0, max = 1.0)    
    # params.add('offset', value = p0[3], min = 0.0, max = 1.0)
    params.add('po', value = p0[0])
    params.add('sigma', value = p0[1])
    params.add('amp', value = p0[2])    
    params.add('offset', value = p0[3])
    out = minimize(Residual, params, args=(x, y, ySEM))
    # out = minimize(Residual, params, args=(x, y, np.ones((len(x), ))))
    return out
    
def PeriodicGaussian(theta, po, sigma, a, offset):
    # po, sigma in degrees
    out = 0.0
    po = AngleInRange(po)
    sigma = AngleInRange(sigma)
    for n in np.arange(-10, 10):
	out += a * np.exp(-(theta - po + 180.0 * n)**2 / (sigma**2)) + offset
    return out

def AngleInRange(theta):
    while theta > 180:
        theta -= 180
    while theta < 0:
        theta += 180
    return theta


def GetPhase(firingRate, atTheta, IF_IN_RANGE = False):
    out = np.nan
    zk = np.dot(firingRate, np.exp(2j * atTheta * np.pi / 180))
    out = np.angle(zk) * 180.0 / np.pi
    if IF_IN_RANGE:
	if(out < 0):
	    out += 360
    return out * 0.5

def StoreToTxt(x, y):
    x.shape = len(x), 1
    y.shape = len(x), 1    
    np.savetxt('data.txt', np.concatenate((x, y), axis = 1))



def FitPeriodicGaussian(x, y, ySEM, p0 = [90, 30, 0.1, 1], IF_PLOT_FIT = False, MIN_PEAK_RATE = 0.0):
    # x = np.concatenate((x, [np.pi]))
    # y = np.concatenate((y, [y[0]]))
    # ySEM = np.concatenate((ySEM, [ySEM[0]]))
    fitParams = np.empty((4, ))
    fitParams[:] = np.nan
    fitError = np.nan
    chiSquare = np.nan
    IS_RESPONSIVE = False
    if(np.max(y) > MIN_PEAK_RATE):
	IS_RESPONSIVE = True
	nPhis = len(x)
	nParams = 4
	degsOfFreedom = nPhis - nParams
	significanceVal = 0.05
	po = GetPhase(y, x * np.pi / 180.0, IF_IN_RANGE = True)
	# p0[0] = x[np.argmax(y)]
	# print po * 180 / np.pi
	# p0[0] = po
	# p0[2] = np.nanmax(y)
	# plt.vlines(po * 180. / np.pi, 0, np.nanmax(y))

	# print 'x:', x
	# print 'y:', y
	try:
    	    theta = np.linspace(0, np.pi, 100)
	    bnds = ((-np.inf, -np.inf, 0, 0), (np.inf, np.inf, np.inf, 1))
	    fitParams, fitError = curve_fit(PeriodicGaussian, x, y, p0 = p0, max_nfev = 4000, bounds = bnds)	    
	    fitY = PeriodicGaussian(x, *fitParams)
	    vidx = ySEM > 0
	    # print fitParams
   	    # print AngleInRange(fitParams[0]), AngleInRange(fitParams[1])    
            # chiSquare = ((y[vidx] - fitY[vidx])**2 / ySEM[vidx]**2).sum()
            # chiSquare = (((y[vidx] - fitY[vidx]) / ySEM[vidx])**2).sum()
	    chiSquare = ((y[vidx] - fitY[vidx])**2 / ySEM[vidx]).sum()
    	    # chiSquare = ((y[vidx] - fitY[vidx])**2 / fitY[vidx]).sum()
	    #################################################################
	    # probability the chi^2 is exceedes the computed value just by chance
	    qVal = GammaQ(degsOfFreedom, chiSquare)
	    # print qVal
	    if qVal < 1:
		fitParamsList = []
		qvalList = []
		chi2List = []
		tmpErrorList = []
		for kk in range(10):
		    # p0 = [x[np.argmax(y)], 100 + 50 * np.random.randn(), np.nanmax(y), np.random.rand()]
   		    # p0 = [x[np.argmax(y)], 100 + 50 * np.random.randn(), np.nanmax(y) + np.random.rand(), 1.0]
		    p0 = [x[np.argmax(y)], 50 + 20 * np.random.randn(), np.nanmax(y) + np.random.rand(), 1.0]
		    # print 'p0 = ', p0
		    fitParamsTmp, fitErrorTmp = curve_fit(PeriodicGaussian, x, y, p0 = p0, max_nfev = 4000, bounds = bnds)

		    # print '     ', fitParamsTmp
		    fitYTmp = PeriodicGaussian(x, *fitParamsTmp)
		    vidx = ySEM > 0
		    tmpErrorList.append(((y[vidx] - fitYTmp[vidx])**2).sum())
		    # chiSquareTmp = (((y[vidx] - fitYTmp[vidx]) / ySEM[vidx])**2).sum()
   		    chiSquareTmp = ((y[vidx] - fitYTmp[vidx])**2 / ySEM[vidx]).sum()
		    # chiSquareTmp = ((y[vidx] - fitYTmp[vidx])**2 / fitY[vidx]).sum()
		    chi2List.append(chiSquareTmp)
		    qvalList.append(GammaQ(degsOfFreedom, chiSquareTmp))
		    fitParamsList.append(fitParamsTmp)
   	    	    # print AngleInRange(fitParamsTmp[0]), AngleInRange(fitParamsTmp[1]), chiSquareTmp
	            # plt.plot(theta * 180 / np.pi, PeriodicGaussian(theta * 180 / np.pi, *fitParamsTmp), label = '%s, %.6s'%(kk, chiSquareTmp))
		    # print kk, chiSquareTmp
		# print len(tmpErrorList), len(fitParamsList), np.argmin(tmpErrorList)
		fitParams = fitParamsList[np.argmin(tmpErrorList)]
		chiSquare = chi2List[np.argmin(tmpErrorList)]
		qVal = qvalList[np.argmin(tmpErrorList)]
		# np.argmin(tmpErrorList)
		# plt.legend()
		# print np.argmin(tmpErrorList)
        	
	    # print 'chi = ', chiSquare
	    IS_GOOD_FIT = qVal > significanceVal
	    fitParams[0] = AngleInRange(fitParams[0])
	    fitParams[1] = AngleInRange(fitParams[1])
	    
 	    if IF_PLOT_FIT:
		# print 'xx'*50
		plt.ion()
		# plt.plot(theta * 180 / np.pi, PeriodicGaussian(theta * 180 / np.pi, *fitParams), 'r')
		plt.plot(theta * 180 / np.pi, PeriodicGaussian(theta * 180 / np.pi, *fitParams))
		plt.legend()
		plt.show()
	except RuntimeError, e:
	    fitParams = np.empty((4, ))
	    fitParams[:] = np.nan
	    fitError = np.nan
	    chiSquare = np.nan
	    print e
    # ipdb.set_trace()
    if(~np.any(np.isnan(fitParams))):
	return fitParams, fitError, chiSquare, IS_GOOD_FIT, IS_RESPONSIVE 
        print fitParams
    else:
	return fitParams, fitError, chiSquare, False, False

def ComputeTuningSEM(p, gamma, nPhis, mExt, mExtOne, nChnks = 2, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_PLOT = False, IF_FIT = False, nNeurons = 10000, neuronType = 'E'):
    # plot tc with sem given atleast 2 chunks of simulations
    NE = N
    NI = N
    nNeurons = NE + NI
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    tc = np.empty((nChnks, nNeurons, nPhis))
    tc[:] = np.nan
    for kChunk in range(nChnks):
	for i, iPhi in enumerate(phis):
	    fr = LoadFrChnk(p, gamma, iPhi, mExt, mExtOne, kChunk, trNo, T, NE, K, nPop)
	    tc[kChunk, :, i] = fr
    tcSem = np.squeeze(np.nanstd(tc, 0)) / np.sqrt(nChnks)
    # tcSem = np.squeeze(np.nanstd(tc, 0))
    tc = np.squeeze(np.nanmean(tc, 0))
    IS_GOOD_FIT = np.zeros((nNeurons, ))
    IS_SINGLE_PEAKED = np.zeros((nNeurons, ))
    IS_RESPONSIVE = np.zeros((nNeurons, ))    
    if IF_FIT and not IF_PLOT:
        chiSquareArray = np.empty((nNeurons, ))
	chiSquareArray[:] = np.nan
	fitParams = []
	fitError = []
	for k in range(nNeurons):
	    # print k
	    yy = tc[k, :]
	    yysem = tcSem[k, :]
	    fitParamsTmp, fitErrorTmp, chiSquareArray[k], IS_GOOD_FIT[k], IS_RESPONSIVE[k] = FitPeriodicGaussian(phis, yy, yysem)
	    fitParams.append(fitParamsTmp)
	    fitError.append(fitErrorTmp)
	print tc.shape, tcSem.shape, chiSquareArray.shape, IS_GOOD_FIT.shape, IS_RESPONSIVE.shape
	np.savez_compressed('./data/fit_p%sg%s_m0%s_mOne%s_nPhis%s'%(int(10 * p), int(10 * gamma), int(1e3 * mExt), int(1e3 * mExtOne), nPhis), tc=tc, tcSem=tcSem, chiSquareArray=chiSquareArray, IS_GOOD_FIT=IS_GOOD_FIT, IS_RESPONSIVE=IS_RESPONSIVE, fitParams = fitParams)
	# np.save('./data/fit_p%sg%s_m0%s_mOne%s'%(int(10 * p), int(10 * gamma), int(1e3 * mExt), int(1e3 * mExtOne)), [tc, tcSem, chiSquareArray, IS_GOOD_FIT, IS_RESPONSIVE])
    if IF_PLOT:
        plt.ion()
    	for i in np.random.randint(0, NE, 101):
	    yy = tc[i, :]
	    yysem = tcSem[i, :]
	    (_, caps, _) = plt.errorbar(phis, yy, fmt = 'ko-', markersize = 3, yerr = yysem, lw = 0.8, elinewidth=0.8)
	    for cap in caps:
		cap.set_markeredgewidth(0.8)
	    osi = sc.OSI(tc[i, :], phis)
	    po = GetPhase(tc[i, :], phis, IF_IN_RANGE = True)
            IS_SINGLE_PEAK = False
	    print yy
	    if IF_FIT:
		fitParams, _, chiSquare, IS_SINGLE_PEAK, _ = FitPeriodicGaussian(phis, yy, yysem, IF_PLOT_FIT = IF_PLOT)
         	print 'neuron#', i, 'firparams = ', fitParams
		plt.show()
		if(~np.any(np.isnan(fitParams))):
		    spVal = GammaQ(nPhis - 4, chiSquare) > 0.75
		    plt.title(r'$neuron# %s $'%(i, osi))		    
                    # plt.title(r'$\mathrm{neuron} %s, \, osi = %.4s, \, \chi^2 = %.4s, \, sp=%s$'%(i, osi, chiSquare, spVal))
        	else:
		    plt.title(r'$neuron# %s, \, osi = %.4s, \,$'%(i, osi))
	    # plt.ion()
	    # plt.show()
	    # plt.waitforbuttonpress()
	    plt.clf()

    if IF_FIT and not IF_PLOT:
	return tc, tcSem, chiSquareArray, IS_GOOD_FIT, IS_RESPONSIVE
    else:
	return tc, tcSem


def PrintTuningBook(p, gamma, mExt, mExtOne, nPhis,  neuronType, nNeurons, fname, trNo = 0, nChnks= 2, T = 1000, NE = 10000, NI = 10000, K = 1000, nPop = 2, IF_FIT = True, IF_PLOT = True, IF_VERBOSE = True, IF_GEN_RAND_NEURONS = False):
    print p, gamma, mExt, mExtOne, nPhis,  neuronType, nNeurons, fname
    doc = Document(fname)
    # doc.packages.append(Package('geometry', options=['left=2cm', 'right=2cm']))
    doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=.8cm', 'bottom=1.2cm']))
    width = r'.3\linewidth'
    nFigsPerPage = 18
    print nNeurons, nFigsPerPage
    nPages = int(np.ceil(nNeurons / float(nFigsPerPage)))
    plt.figure()
    plt.ioff()
    theta = np.linspace(0, 180, nPhis, endpoint = False)
    # doc = Document(fname)
    if IF_GEN_RAND_NEURONS:
	if neuronType == 'E':
	    rndNeurons = np.random.randint(0, NE, nNeurons)
	else:
	    rndNeurons = np.random.randint(NE, NE + NI, nNeurons)
        np.save('./data/twopop/neuronIdx_n100_NE', rndNeurons)
    else:
	rndNeurons = np.load('./data/twopop/neuronIdx_n100_NE.npy')

    rndNeurons = np.sort(rndNeurons)
    nNeurons = NE + NI
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    tc = np.empty((nChnks, nNeurons, nPhis))
    tc[:] = np.nan
    
    # tcLast = np.empty((nNeurons, nPhis))
    # tcLast[:] = np.nan
    # for i, iPhi in enumerate(phis):
    # 	fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, trNo, T, NE, K, nPop)
    # 	tcLast[:, i] = fr
    
    for kChunk in range(nChnks):
	for i, iPhi in enumerate(phis):
	    fr = LoadFrChnk(p, gamma, iPhi, mExt, mExtOne, kChunk, trNo, T, NE, K, nPop)
	    tc[kChunk, :, i] = fr
    tcSem = np.squeeze(np.nanstd(tc, 0)) / np.sqrt(nChnks)
    # tcSem = np.squeeze(np.nanstd(tc, 0))
    tc = np.squeeze(np.nanmean(tc, 0))
    fitError = []    
    for kk in range(nPages):
	rndNeuronsx = rndNeurons[kk * nFigsPerPage: (kk + 1) * nFigsPerPage]
	print kk * nFigsPerPage, (kk + 1) * nFigsPerPage
	with doc.create(Figure(position='htbp')) as plot:
	    for i in rndNeuronsx:
		# i = 7304
		yy = tc[i, :]
		yysem = tcSem[i, :]
		# if(np.argmax(yy) == 0):
		#     yy = np.roll(yy, 4)
		#     yysem = np.roll(tcSem[i, :], 4)
		(_, caps, _) = plt.errorbar(phis, yy, fmt = 'ko-', markersize = 3, yerr = yysem, lw = 0.8, elinewidth=0.8)
		# plt.plot(phis, tcLast[i, :], 'c')
		for cap in caps:
		    cap.set_markeredgewidth(0.8)
		osi = sc.OSI(tc[i, :], phis)
		po = GetPhase(tc[i, :], phis, IF_IN_RANGE = True)
		IS_SINGLE_PEAK = False
		if IF_FIT:
		    fitParams, _, chiSquare, IS_SINGLE_PEAK, IS_RESPONSIVE = FitPeriodicGaussian(phis, yy, yysem, IF_PLOT_FIT = IF_PLOT, MIN_PEAK_RATE = 0.)
		    fity = PeriodicGaussian(np.linspace(0, 180, nPhis, endpoint = False), *fitParams)
		    fitMax = np.max(PeriodicGaussian(np.linspace(0, 180, 100, endpoint = 0), *fitParams))
		    tcMax = np.max(yy)
		    fitError = np.sqrt(np.mean((tc[i, :] - fity) **2)) / fitMax
		    peakError = 100.0 * np.abs(tcMax - fitMax) / tcMax
		if IF_VERBOSE:
		    print 'neuron#', i, 'firparams = ', fitParams, 'qprob = ', GammaQ(nPhis - 4, chiSquare)
		# plt.show()
		if(~np.any(np.isnan(fitParams))):
                    # plt.title(r'$\mathrm{neuron:} %s, \, dist = %.6s$'%(i, fitError))
                    plt.title(r'$\mathrm{neuron:} %s, \,\chi^2 =  %s$'%(i, chiSquare))		    
    		    # spVal = GammaQ(nPhis - 4, chiSquare) > 0.05
                    # plt.title(r'$\mathrm{neuron:} %s, \, osi = %.4s, \, \chi^2 = %.6s\, SP = %s$'%(i, osi, chiSquare, int(spVal)))
                    # plt.title(r'$\mathrm{neuron:} %s, \, \chi^2 = %.6s$'%(i, chiSquare))
		    # plt.title(r'$\mathrm{neuron:} %s, \, osi = %.4s, \, \chi^2 = %.6s$'%(i, osi, chiSquare))
                    # plt.title(r'$\mathrm{neuron} %s, \, osi = %.4s$'%(i, osi))
        	else:
                    # plt.title(r'$\mathrm{neuron:} %s, \, osi = %.4s, \,$'%(i, osi))
                    plt.title(r'$\mathrm{neuron:} %s$'%(i,))
		# plt.draw()
		# ipdb.set_trace()
		# pplot = np.max(yy) > 0.18
		if True: #fitError < 0.06 and peakError < 5:
		    with doc.create(SubFigure(position='b', width=NoEscape(width))) as figure:
			figure.add_plot(width=NoEscape(r'\linewidth'), dpi = 300) #*args, **kwargs)
		plt.ion()
	        plt.waitforbuttonpress()
		plt.clf()
    doc.generate_pdf(clean_tex=False)



def GetSinglePeakNeurons(p, gamma, mExt, mExtOne, nPhis,  neuronType, nNeurons, fname, trNo = 0, nChnks= 2, T = 1000, NE = 10000, NI = 10000, K = 1000, nPop = 2, IF_FIT = True, IF_PLOT = True, IF_VERBOSE = True, IF_GEN_RAND_NEURONS = False):
    theta = np.linspace(0, 180, nPhis, endpoint = False)
    if IF_GEN_RAND_NEURONS:
	if neuronType == 'E':
	    rndNeurons = np.random.randint(0, NE, np.min([10 * nNeurons, NE]))
	else:
	    rndNeurons = np.random.randint(NE, NE + NI, np.min([10 * nNeurons, NI]))
        np.save('./data/twopop/neuronIdx_n100_NE', rndNeurons)
    else:
	rndNeurons = np.load('./data/twopop/neuronIdx_n100_NE.npy')
    rndNeurons = np.sort(rndNeurons)
    nNeurons = NE + NI
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    tc = np.empty((nChnks, nNeurons, nPhis))
    tc[:] = np.nan
    for kChunk in range(nChnks):
	for i, iPhi in enumerate(phis):
	    fr = LoadFrChnk(p, gamma, iPhi, mExt, mExtOne, kChunk, trNo, T, NE, K, nPop)
	    tc[kChunk, :, i] = fr
    tcSem = np.squeeze(np.nanstd(tc, 0)) / np.sqrt(nChnks)
    tc = np.squeeze(np.nanmean(tc, 0))
    fitError = []
    singlePeakedIdx = []
    neuronCount = 0
    iterCnt = 0
    print '\/'*10
    print 'fitting params ...'
    while(neuronCount < nNeurons and iterCnt < len(rndNeurons)):
    	i = rndNeurons[iterCnt]
    	yy = tc[i, :]
    	yysem = tcSem[i, :]
    	po = GetPhase(tc[i, :], phis, IF_IN_RANGE = True)
    	IS_SINGLE_PEAK = False
    	if IF_FIT:
    	    fitParams, _, chiSquare, IS_SINGLE_PEAK, IS_RESPONSIVE = FitPeriodicGaussian(phis, yy, yysem, IF_PLOT_FIT = IF_PLOT, MIN_PEAK_RATE = 0)
    	    fity = PeriodicGaussian(np.linspace(0, 180, nPhis, endpoint = False), *fitParams)
    	    fitMax = np.max(PeriodicGaussian(np.linspace(0, 180, 100, endpoint = 0), *fitParams))
    	    tcMax = np.max(yy)
    	    fitError = np.sqrt(np.mean((tc[i, :] - fity) **2)) / fitMax
    	    peakError = 100.0 * np.abs(tcMax - fitMax) / tcMax
    	    if fitError < 0.06 and peakError < 6:
    		neuronCount += 1
    		singlePeakedIdx.append(i)
    	iterCnt += 1
    # 	print iterCnt, neuronCount
	         
    ### print to pdf ###
    print '\/'*10
    print 'printing pdf...'
    print 'n sp neurons = ', len(singlePeakedIdx)
    # singlePeakedIdx = range(2)
    fname = './figs/twopop/p%sg%s_m%s_nphis%s'%(int(p*10), int(gamma*10), int(mExtOne * 1e3), nPhis)
    if len(singlePeakedIdx) > 0:
	doc = Document(fname)
	doc.packages.append(Package('geometry', options=['left=1cm', 'right=1cm', 'top=.8cm', 'bottom=1.2cm']))
	width = r'.3\linewidth'
	nFigsPerPage = 18
	nNeurons = len(singlePeakedIdx)
	nPages = int(np.ceil(nNeurons / float(nFigsPerPage)))
	print 'nPasges = ', nPages
	# plt.ion()
	plt.clf()
	# plt.show()
        # ipdb.set_trace()		
	for kk in range(nPages):
	    singlePeakedIdxList = singlePeakedIdx[kk * nFigsPerPage: (kk + 1) * nFigsPerPage]
	    with doc.create(Figure(position='htbp')) as plot:
		for i in singlePeakedIdxList:
		    # plt.ion()
		    # plt.waitforbuttonpress()
		    
		    yy = tc[i, :]
		    yysem = tcSem[i, :]
		    (_, caps, _) = plt.errorbar(phis, yy, fmt = 'ko-', markersize = 3, yerr = yysem, lw = 0.8, elinewidth=0.8)
		    for cap in caps:
			cap.set_markeredgewidth(0.8)
		    IS_SINGLE_PEAK = False
		    if IF_FIT:
			fitParams, _, chiSquare, IS_SINGLE_PEAK, IS_RESPONSIVE = FitPeriodicGaussian(phis, yy, yysem, IF_PLOT_FIT = IF_PLOT, MIN_PEAK_RATE = 0.)
		    if IF_VERBOSE:
			print 'neuron#', i, 'firparams = ', fitParams, 'qprob = ', GammaQ(nPhis - 4, chiSquare)
		    if(~np.any(np.isnan(fitParams))):
			# plt.title(r'$\mathrm{neuron:} %s, \, \chi^2 = %.4s$'%(i, chiSquare))
			plt.title(r'$\mathrm{neuron:} %s$'%(i, ))			
		    else:
			plt.title(r'$\mathrm{neuron:} %s$'%(i,))
		    with doc.create(SubFigure(position='b', width=NoEscape(width))) as figure:
			    figure.add_plot(width=NoEscape(r'\linewidth'), dpi = 300) #*args, **kwargs)
		    # plt.ion()
		    # plt.waitforbuttonpress()
		    plt.clf()
# ;        ipdb.set_trace()	
	doc.generate_pdf(clean_tex=False)


                         

    
    
def PrintTuningBookManual(p, gamma, mExt, mExtOne, nPhis,  neuronType, nNeurons, fname, trNo = 0, nChnks= 2, T = 1000, NE = 10000, NI = 10000, K = 1000, nPop = 2, IF_FIT = True, IF_PLOT = True, IF_VERBOSE = True):
    print p, gamma, mExt, mExtOne, nPhis,  neuronType, nNeurons, fname
    theta = np.linspace(0, 180, nPhis, endpoint = False)
    if neuronType == 'E':
	rndNeurons = np.random.randint(0, NE, nNeurons)
    else:
	rndNeurons = np.random.randint(NE, NE + NI, nNeurons)
    nNeurons = NE + NI
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    tc = np.empty((nChnks, nNeurons, nPhis))
    tc[:] = np.nan
    for kChunk in range(nChnks):
	for i, iPhi in enumerate(phis):
	    fr = LoadFrChnk(p, gamma, iPhi, mExt, mExtOne, kChunk, trNo, T, NE, K, nPop)
	    tc[kChunk, :, i] = fr
    tcSem = np.squeeze(np.nanstd(tc, 0)) / np.sqrt(nChnks)
    tc = np.squeeze(np.nanmean(tc, 0))
    rndNeuronsx = rndNeurons[kk * nFigsPerPage: (kk + 1) * nFigsPerPage]
    for i in rndNeuronsx:
	FIT_DONE = False
	yy = tc[i, :]
	yysem = tcSem[i, :]
	print "neuron# ", i
	while(FIT_DONE == 'n' or FIT_DONE == 0):
	    (_, caps, _) = plt.errorbar(phis, yy, fmt = 'ko-', markersize = 3, yerr = yysem, lw = 0.8, elinewidth=0.8)
	    plt.title(r'$\mathrm{neuron:} %s$'%(i,))
	    po = GetPhase(yy, phis * np.pi / 180.0, IF_IN_RANGE = True)
	    if nIters == 0:
		p0 = [90, 80, 0.1, 1]
		p0[0] = po
		p0[2] = np.nanmax(yy)
	    else:
		fitPO = input('po = ')
		fitSigma = input('sigma = ')
		fitAmp = input('amplitude = ')
		fitOffset = input('offset = ')
		p0 = [fitPO, fitSigma, fitAmp, fitOffset]
		print 'p0 = ', p0
	    for cap in caps:
		cap.set_markeredgewidth(0.8)
	    osi = sc.OSI(tc[i, :], phis)
	    po = GetPhase(tc[i, :], phis, IF_IN_RANGE = True)
	    IS_SINGLE_PEAK = False
	    if IF_FIT:
		fitParams, _, chiSquare, IS_SINGLE_PEAK, IS_RESPONSIVE = FitPeriodicGaussian(phis, yy, yysem, IF_PLOT_FIT = IF_PLOT, MIN_PEAK_RATE = 0, p0 = p0)
	    if IF_VERBOSE:
		print 'neuron#', i, 'firparams = ', fitParams, 'qprob = ', GammaQ(nPhis - 4, chiSquare)
	    if(~np.any(np.isnan(fitParams))):
		plt.title(r'$\mathrm{neuron:} %s, \, \chi^2 = %.6s$'%(i, chiSquare))
	    else:
		plt.title(r'$\mathrm{neuron:} %s$'%(i,))
	    pplot = np.max(yy) > 0.18
	    FIT_DONE = raw_input("Accept fit? ")
	    print FIT_DONE
	    nIters += 1
	    if FIT_DONE == 1 or FIT_DONE == 'y':
		plt.clf()
		(_, caps, _) = plt.errorbar(phis, yy, fmt = 'ko-', markersize = 3, yerr = yysem, lw = 0.8, elinewidth=0.8)
		fitParams, _, chiSquare, IS_SINGLE_PEAK, IS_RESPONSIVE = FitPeriodicGaussian(phis, yy, yysem, IF_PLOT_FIT = IF_PLOT, MIN_PEAK_RATE = 0, p0 = p0)
		theta = np.linspace(0, 180, 100)
		plt.plot(theta, PeriodicGaussian(theta, *fitParams), 'r')
		with doc.create(SubFigure(position='b', width=NoEscape(width))) as figure:
		    figure.add_plot(width=NoEscape(r'\linewidth'), dpi = 300) #*args, **kwargs)
		    plt.clf()
    		    # ipdb.set_trace()
		plt.clf()
    doc.generate_pdf(clean_tex=False)

if __name__ == "__main__":
    p = float(sys.argv[1])
    gamma = float(sys.argv[2])
    mExt = float(sys.argv[3])
    mExtOne = float(sys.argv[4])
    IF_PRINT_BOOK_OF_TUNING = int(sys.argv[5])
    nPhis = 16
    trialNo = 10
    
    ComputeTuningSEM(p, gamma, nPhis, mExt, mExtOne, IF_PLOT = 0, IF_FIT = 1, trNo = trialNo)
    
