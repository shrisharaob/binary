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

def LoadFr(p, gamma, phi, mExt, mExtOne, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False):
    baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'

    if nPop == 1:
    	baseFldr = baseFldr + 'onepop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if nPop == 2:
    	baseFldr = baseFldr + 'twopop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if IF_VERBOSE:
    	print baseFldr

    # if nPop == 1:
    # 	baseFldr = baseFldr + 'onepop/data/old/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    # if nPop == 2:
    # 	baseFldr = baseFldr + 'twopop/data/old/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    # if IF_VERBOSE:
    # 	print baseFldr
	
    filename = 'meanrates_theta%.6f_tr%s.txt'%(phi, trNo)
    return np.loadtxt(baseFldr + filename)

def LoadFrChnk(p, gamma, phi, mExt, mExtOne, chkId, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False):
    baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'

    if nPop == 1:
    	baseFldr = baseFldr + 'onepop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if nPop == 2:
    	baseFldr = baseFldr + 'twopop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if IF_VERBOSE:
    	print baseFldr 

    # if nPop == 1:
    # 	baseFldr = baseFldr + 'onepop/data/old/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    # if nPop == 2:
    # 	baseFldr = baseFldr + 'twopop/data/old/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    # if IF_VERBOSE:
    # 	print baseFldr
	
    filename = 'meanrates_theta%.6f_tr%s_chnk%s.txt'%(phi, trNo, chkId)
    return np.loadtxt(baseFldr + filename)

 
def LoadSpkTimes(p, gamma, phi, mExt, mExtOne, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2):
    baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
    if nPop == 1:
	baseFldr = baseFldr + 'onepop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if nPop == 2:
	baseFldr = baseFldr + 'twopop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    print 'loading from fldr: ', baseFldr
    return np.loadtxt(baseFldr + 'spktimes_theta%.6f_tr%s.txt'%(phi, trNo), delimiter = ';')

def SpikesInInterval(st, spkStart, spkEnd = -1):
    if(spkEnd == -1):
        spkEnd = np.max(st[:, 1])
        idx = st[:, 1] >= spkStart
    else:
        idx = np.logical_and(st[:, 1] >= spkStart, st[:, 1] <= spkEnd)
    return st[idx, :]

def PlotRasterBasic(st):
    plt.ion()
    plt.plot(st[:, 1], st[:, 0], '|k')
    plt.draw()

def PlotRaster(st, spkStart, spkStop, NE, NI, IF_COLOR = False):
    plt.ion()
    st = SpikesInInterval(st, spkStart, spkStop)
    if IF_COLOR:
        idxNE = st[:, 0] < NE
        idxNI = st[:, 0] > NE
        stNE = st[idxNE, :]
        stNI = st[idxNI, :]
        plt.plot(stNE[:, 1], stNE[:, 0], '|k')
        plt.plot(stNI[:, 1], stNI[:, 0], '|', color = [0.7, 0.1, 0.1])
        plt.draw()
        print 'hello'
    else:
        plt.plot(st[:, 1], st[:, 0], '.k')
        plt.draw()
    plt.xlabel('Time')
    plt.ylabel('Neuron #')

def RasterPlot(st, spkStart = 0, spkStop = -1, ):
    st = SpikesInInterval(st, spkStart, spkStop)
    print st.shape
    neuronIdx = np.unique(st[:, 0])
    nNeurons = np.size(neuronIdx)
    nSpks, _ = st.shape
    totalvLength = 200.0
    vLenght = 100.0
    vLineOffset = (totalvLength - vLenght) * 0.5
    plt.figure()
    plt.ion()
    x = np.array([])
    y = np.array([])
    for idx, iNeuron in enumerate(neuronIdx):
        iSpkTimes = st[st[:, 1] == iNeuron, 0]
        x = np.r_[x, iSpkTimes]
        y = np.r_[y, totalvLength * iNeuron * np.ones((np.size(iSpkTimes),)) + vLineOffset - vLenght]
    plt.vlines(x, y, y + vLenght)
    plt.ylim(y[0] - 1 - vLenght, y[-1] + 1 + vLenght)
    plt.yticks(neuronIdx * totalvLength)
    plt.gca().set_yticklabels(neuronIdx)
    plt.draw()
    plt.show()

def AutoCorr(x, corrLength = "same"):
    # x : spike Times in ms
    N = len(x)
    nPointFFT = int(np.power(2, np.ceil(np.log2(len(x)))))
    fftX = np.fft.fft(x, nPointFFT)
    return np.abs(np.abs(np.fft.ifft(np.multiply(fftX, np.conj(fftX)))))

def AvgAutoCorrInInterval(starray, neuronsList, spkTimeStart, spkTimeEnd, minSpks = 100, maxTimeLag = 100, NE = 10000, NI = 10000, NFF = 10000, fileTag = 'E', theta = 0, TAU_FF = 4, TAU_E = 4, TAU_I = 2):
    N = len(neuronsList)
    simDT = (TAU_FF * TAU_E * TAU_I) / (NE * TAU_FF * TAU_I + NI * TAU_FF * TAU_E + NFF * TAU_E * TAU_I)
    pcDone = 0
    simDuration = spkTimeEnd - spkTimeStart
    nTimeLagBins = int(2 * maxTimeLag) # works if downsample bin size = 1ms otherwise divide by value
    avgCorr = np.zeros((int(np.power(2, np.ceil(np.log2(simDuration + simDT)))), ))
    binFlag = 0;
    nValidNeurons = 0;
    downSampleBinSize = 1
    spkBins = np.arange(spkTimeStart, spkTimeEnd + simDT, downSampleBinSize)
    nSpkBins = len(spkBins)
    avgRate = 0
    popMeanRate = [];
    for i, kNeuron in enumerate(neuronsList):
        spkTimes = starray[starray[:, 0] == kNeuron, 1]
        spksTimes = spkTimes[spkTimes > spkTimeStart]
        meanRate = float(spkTimes.size) / float(simDuration)
        nValidNeurons += 1
        avgRate += meanRate
        # print spkTimes.size, spkTimes.size >= minSpks
        if(spkTimes.size >= minSpks):
            # print meanRate
            st = np.histogram(np.squeeze(spkTimes), spkBins)
            tmpCorr = AutoCorr(st[0])
            avgCorr += tmpCorr / ((downSampleBinSize) **2 * nSpkBins * meanRate)
    avgCorr = avgCorr / nValidNeurons
    bins = np.array(downSampleBinSize)
    avgCorr[np.argmax(avgCorr)] = 0.0
    avgCorr = avgCorr[: maxTimeLag]
    print 'avgRate = ', avgRate / N
    if(len(avgCorr) > 0):
        return avgCorr
    else :
        return 0

def PlotAC(st, nNeurons, spkTimeStart, spkTimeEnd, minSpks = 100, maxTimeLag = 100, NE = 10000, NI = 10000, NFF = 10000, fileTag = 'E', theta = 0, TAU_FF = 4, TAU_E = 4, TAU_I = 2):
    nENeurons = np.arange(NE)
    nINeurons = np.arange(NE, NE + NI)
    if nNeurons < NE:
        nENeurons = np.random.randint(0, NE, nNeurons)
    if nNeurons < NI:
        nINeurons = np.random.randint(NE, NE + NI, nNeurons)
    acE = AvgAutoCorrInInterval(st, nENeurons, spkTimeStart, spkTimeEnd, NE = NE, NI = NI, NFF = NFF, TAU_E = TAU_E, TAU_I = TAU_I, TAU_FF = TAU_FF, theta = theta, maxTimeLag = maxTimeLag)
    acI = AvgAutoCorrInInterval(st, nINeurons, spkTimeStart, spkTimeEnd, NE = NE, NI = NI, NFF = NFF, TAU_E = TAU_E, TAU_I = TAU_I, TAU_FF = TAU_FF, theta = theta, maxTimeLag = maxTimeLag, minSpks = minSpks)
    plt.figure()
    plt.plot(acE, 'k', label = 'E')
    plt.plot(acI, color = [0.7, 0.1, 0.1], label = 'I')
    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 10})
    plt.xlabel('Time')
    # plt.ylabel('')
    plt.title('AC')
    return acE, acI
    
def CV(spkTimes):
    out = np.array([np.nan, np.nan])
    if(spkTimes.size > 3):
        isi = np.diff(spkTimes)
        cv = isi.std() / isi.mean()
        cv2 = CV2(isi)
        out = np.array([cv, cv2])
    return out

def CV2(isi):
    denom = isi[1:] + isi[:-1]
    cv2 = 2.0 * np.abs(np.diff(isi)) / denom
    return np.nanmean(cv2)

def CVofList(starray, listOfneurons, spkTimeStart = 200):
    out = np.empty((len(listOfneurons), 2))
    out[:] = np.nan
    for i, iNeuron in enumerate(listOfneurons):
        spkTimes = starray[starray[:, 0] == iNeuron, 1]
        spksTimes = spkTimes[spkTimes > spkTimeStart]
        out[i, :] = CV(spkTimes)
    return out #out[:, 0] is CV, out[:, 1] is CV2

def PlotCVDistr(starray, NE, NI, nNeurons, spkTimeStart = 500, nBins = 10):
    plt.figure()
    print 'computeing CV, E ... ',
    sys.stdout.flush()
    cvE = CVofList(starray, np.random.randint(0, NE, nNeurons), spkTimeStart)
    print 'done'
    print 'computeing CV, I ... ',    
    cvI = CVofList(starray, np.random.randint(NE, NE + NI, nNeurons), spkTimeStart)
    print 'done'
    avgCVE = np.nanmean(cvE[:, 0])
    avgCVI = np.nanmean(cvI[:, 0])
    avgCV2E = np.nanmean(cvE[:, 1])
    avgCV2I = np.nanmean(cvI[:, 1])
    print 'avg CVE = ', avgCVE, 'avg CVI = ', avgCVI
    print 'avg CV2E = ', avgCV2E, 'avg CV2I = ', avgCV2I    
    plt.hist(cvE[~np.isnan(cvE[:, 0]), 0], nBins, histtype = 'step', normed = 1, label = 'E', color = 'k')
    plt.hist(cvI[~np.isnan(cvI[:, 0]), 0], nBins, histtype = 'step', normed = 1, label = 'I', color = [0.7, 0.1, 0.1])
    _, ymax = plt.ylim()
    plt.vlines(avgCVE, 0, ymax, color = 'k', lw = 2)
    plt.vlines(avgCVI, 0, ymax, color = [0.7, 0.1, .1], lw = 2)    
    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 10})
    plt.xlabel('CV')
    plt.ylabel('Density')
    plt.figure()
    plt.hist(cvE[~np.isnan(cvE[:, 1]), 1], nBins, histtype = 'step', normed = 1, label = 'E', color = 'k')
    plt.hist(cvI[~np.isnan(cvI[:, 1]), 1], nBins, histtype = 'step', normed = 1, label = 'I', color = [0.7, 0.1, 0.1])
    _, ymax = plt.ylim()    
    plt.vlines(avgCV2E, 0, ymax, color = 'k', lw = 2)
    plt.vlines(avgCV2I, 0, ymax, color = [0.7, 0.1, .1], lw = 2)    
    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 10})
    plt.xlabel(r'$CV_{2}$')
    plt.ylabel('Density')
        
def PlotInstantRates(stArray, spkTimeStart, spkTimeEnd, NE, NI, NFF, windowSize, IF_NEW_FIG = True, TAU_FF = 4.0, TAU_E = 4.0, TAU_I = 2.0):
    st =  SpikesInInterval(stArray, spkTimeStart, spkTimeEnd)
    dt = (TAU_FF * TAU_E * TAU_I) / (NE * TAU_FF * TAU_I + NI * TAU_FF * TAU_E + NFF * TAU_E * TAU_I)
    windowSize *= dt
    print dt
    binEdges = np.arange(spkTimeStart, spkTimeEnd + windowSize, windowSize)
    if IF_NEW_FIG:
        plt.figure()
    spkTimesE = st[st[:, 0] < NE, 1]
    spkTimesI = st[st[:, 0] >= NE, 1]
    cntE, binsE = np.histogram(spkTimesE, binEdges)
    cntI, binsI = np.histogram(spkTimesI, binEdges)
    print binsE.shape, binsE[-1], binsE[0]
    print cntE
    print cntI
    plt.plot(binsE[:-1], cntE / (float(windowSize) * NE), color = 'k')
    plt.plot(binsI[:-1], cntI / (float(windowSize) * NI), color = [0.7, .1, .1])
    plt.grid()
    
def OSI(firingRate, atTheta):
    out = np.nan
    validIdx = firingRate > 0
    # zk = 2.0 * np.dot(firingRate[validIdx], np.exp(2j * atTheta[validIdx] * np.pi / 180)) / validIdx.sum()
    # print np.absolute(zk)
    mA1 = M1Component(firingRate) * 0.5 
    mA0 = np.nanmean(firingRate)
    # print mA1, mA0
    if(mA0 > 0):
        out = mA1 / mA0
    if(out > 1):
        out = np.nan
    return out

def GetPhase(firingRate, atTheta, IF_IN_RANGE = False):
    out = np.nan
    zk = np.dot(firingRate, np.exp(2j * atTheta * np.pi / 180))
    out = np.angle(zk) * 180.0 / np.pi
    if IF_IN_RANGE:
	if(out < 0):
	    out += 360
    return out * 0.5

def OSIOfPop(firingRates, atThetas):
    # thetas in degrees
    nNeurons, nThetas = firingRates.shape
    out = np.zeros((nNeurons, ))
    for i in range(nNeurons):
        out[i] = OSI(firingRates[i , :], atThetas)
    return out


def POofPopulation(tc, theta = np.arange(0.0, 180.0, 22.5)):
    # return value in degrees
    nNeurons, _ = tc.shape
    po = np.zeros((nNeurons, ))
    for kNeuron in np.arange(nNeurons):
        po[kNeuron] = GetPhase(tc[kNeuron, :], theta)
    return po 


def PltOSIHist(p, gamma, nPhis, mExt, mExtOne, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_NEW_FIG = True, color = 'k'):
    NE = N
    NI = N
    tc = np.zeros((NE + NI, nPhis))
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    if IF_NEW_FIG:
	plt.figure()
    for i, iPhi in enumerate(phis):
        print 'loading from fldr: ', 
	fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, trNo, T, NE, K, nPop, IF_VERBOSE = True)
        if(len(fr) == 1):
            if(np.isnan(fr)):
                print 'file not found!'
	tc[:, i] = fr
    osi = OSIOfPop(tc[:NE, :], phis)
    # plt.xlabel(r"$\mathrm{OSI} \,\,\,\,  (m_{E, i}^{(1)})$")
    plt.xlabel('OSI', fontsize = 12)    
    plt.ylabel('Density', fontsize = 12)
    plt.hist(osi[~np.isnan(osi)], 27, normed = 1, histtype = 'step', label = r'$p = %s,\, \gamma = %s, \, m_0^{(1)} = %s$'%(p, gamma, mExtOne), color = color, lw = 2)
    # plt.title(r'$N = %s,\, K = %s,\, m_0^{(0)} = %s,\, m_0^{(1)} = %s$'%(NE, K, mExt, mExtOne))
    plt.title(r'$N = %s,\, K = %s,\, m_0^{(0)} = %s$'%(NE, K, mExt,))
    _, ymax = plt.ylim()
    plt.vlines(np.nanmean(osi), 0, ymax, lw = 2, color = color)
    print "mean OSI = ", np.nanmean(osi)
    # print "MEAN OSI = ", np.nanmean(osi)
    return osi

def CompareOSIHist(pList, gList, nPhis, mExt, mExtOneList, trNo, N = 10000, K = 1000, nPop = 2, T = 1000, IF_NEW_FIG = True, clrCntr = 0):
    if IF_NEW_FIG:
	plt.figure()
    colors = [plt.cm.Dark2(i) for i in np.linspace(0, 1, 1 + clrCntr + len(pList) * len(gList) * len(mExtOneList), endpoint = False)]
    for mExtOne in mExtOneList:
	for p in pList:
	    for gamma in gList:
		try:
		    PltOSIHist(p, gamma, nPhis, mExt, mExtOne, trNo = trNo, IF_NEW_FIG = False, color = colors[clrCntr])
		    clrCntr += 1
		except IOError:
		    print "p = ", p, " gamma = ", gamma, " trial# ", trNo, " file not found"
    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 10})
    plt.gca().set_position([0.2, 0.2, .65, .65])
    if nPop == 2:
	plt.savefig("./figs/twopop/compareOSI_.png")
    
def ComputeTuningForNeuron(p, gamma, nPhis, mExt, mExtOne, neuronIdx, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000):
    NE = N
    NI = N
    tc = np.zeros((nPhis))
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    for i, iPhi in enumerate(phis):
	fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, trNo, T, NE, K, nPop)
	# print fr.shape
	tc[i] = fr[neuronIdx]
    plt.ioff()
    # plt.figure()
    plt.plot(phis, tc, 'ks-')
    osi = OSI(tc, phis)
    plt.title('neuron# %s, osi = %s'%(neuronIdx, osi))
    return tc, plt.gca()

def ComputeTuningForNeuronSEM(p, gamma, nPhis, mExt, mExtOne, neuronIdx, nChnks = 2, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000):
    NE = N
    NI = N
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    tc = np.empty((2, nPhis))
    tc[:] = np.nan
    tmp = np.empty((nPhis))
    tmp[:] = np.nan
    # ipdb.set_trace()
    for kChunk in range(nChnks):
	for i, iPhi in enumerate(phis):
	    fr = LoadFrChnk(p, gamma, iPhi, mExt, mExtOne, kChunk, trNo, T, NE, K, nPop)
	    tmp[i] = fr[neuronIdx]
	tc[kChunk, :] = tmp
    tc = np.nanmean(tc, 0)
    tcSem = np.nanstd(tc, 0) / np.sqrt(nChnks)
    # plt.ioff()
    (_, caps, _) = plt.errorbar(phis, tc, fmt = 'ko-', markersize = 3, yerr = tcSem, lw = 0.8, elinewidth=0.8)
    for cap in caps:
	# cap.set_color('red')
	cap.set_markeredgewidth(0.8)
    osi = OSI(tc, phis)
    plt.title('neuron# %s, osi = %s'%(neuronIdx, osi))
    return tc, plt.gca()

def ComputeTuningSEM(p, gamma, nPhis, mExt, mExtOne, nChnks = 2, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_PLOT = False, IF_FIT = False, nNeurons = 10000, neuronType = 'E'):
    # plot tc with sem given atleast 2 chunks of simulations
    NE = N
    NI = N
    nNeurons = NE + NI
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    tc = np.empty((2, nNeurons, nPhis))
    tc[:] = np.nan
    nParams = 4
    degsOfFreedom = nPhis - nParams
    significanceVal = 0.05
    criticalChi2 = ChiSquared.isf(significanceVal, degsOfFreedom)
    print 'criticalChi2 = ', criticalChi2
    # ipdb.set_trace()
    for kChunk in range(nChnks):
	for i, iPhi in enumerate(phis):
	    fr = LoadFrChnk(p, gamma, iPhi, mExt, mExtOne, kChunk, trNo, T, NE, K, nPop)
	    tc[kChunk, :, i] = fr
    # ipdb.set_trace()
    tcSem = np.squeeze(np.nanstd(tc, 0)) / np.sqrt(nChnks)
    tcSemOverAngles = np.nanstd(np.nanmean(tc, 2), 0) / np.sqrt(nChnks)
    tc = np.squeeze(np.nanmean(tc, 0))
    # ipdb.set_trace()
    if IF_FIT and not IF_PLOT:
	chiSquareArray = np.empty((nNeurons, ))
	chiSquareArray[:] = np.nan
	avgAbsResidual = np.empty((nNeurons, nPhis))
	avgAbsResidual[:] = np.nan
	for k in range(nNeurons):
	# for k in range(10):
	    _, _, chiSquareArray[k], avgAbsResidual[k, :] = FitVonMisses(phis * np.pi / 180, tc[k, :])
    if IF_PLOT:
        plt.ion()
	
    	for i in np.random.randint(0, NE, 101):
	    # ipdb.set_trace()
	    # i = 2571
	    # i = 1664
	    # plt.plot(np.concatenate((phis, [180])), np.concatenate((tc[i, :], [tc[i, 0]])), 'ks-')
	    (_, caps, _) = plt.errorbar(np.concatenate((phis, [180])), np.concatenate((tc[i, :], [tc[i, 0]])), fmt = 'ko-', markersize = 3, yerr = np.concatenate((tcSem[i, :], [tcSem[i, 0]])), lw = 0.8, elinewidth=0.8)
	    # print tcSem[i, :]
	    for cap in caps:
		# cap.set_color('red')
		cap.set_markeredgewidth(0.8)
	    osi = OSI(tc[i, :], phis)
	    po = GetPhase(tc[i, :], phis, IF_IN_RANGE = True)
	    _, ymax = plt.ylim()
	    plt.vlines(po, 0, ymax, color = 'k')
            IS_SINGLE_PEAK = False
	    if IF_FIT:
                fitParams, _, chiSquare, avgAbsRes = FitVonMisses(phis * np.pi / 180, tc[i, :], IF_PLOT_FIT = True)
		print fitParams

		if(~np.any(np.isnan(fitParams))):
		    # print i, avgAbsRes, chiSquare
		    # viddx = tcSem[i, :] != 0
		    # tmpResidual = avgAbsRes[:-1]
		    # print tmpResidual
		    # print tcSem[i, viddx]
		    # avgAbsRes = np.mean(tmpResidual[viddx] / tcSem[i, viddx])
                    avgAbsRes /= tcSemOverAngles[i]
		    avgAbsRes = avgAbsRes.mean()
		    # print avgAbsRes, chiSquare
		    # fitParams, _, chiSquare, _ = FitVonMisses(phis * np.pi / 180, tc[i, :], IF_PLOT_FIT = True)
		    # ipdb.set_trace()
		    meanRate = np.nanmean(tc[i, :])
		    print 'i=', i, 'avgAbsRes =', avgAbsRes, 'chi2 = ', chiSquare, 'meanR = ', meanRate
		    IS_SINGLE_PEAK = (avgAbsRes <= 2.0) or (chiSquare <= 0.02) and ((np.nanmax(tc[i, :]) - meanRate) >= 0.05)
		    #IF_SINGLE_PEAK(tcSemOverAngles[i], avgAbsRes)
		    # IS_SINGLE_PEAK2 = IF_SINGLE_PEAK(chiSquare, avgAbsRes)
		    if(~np.any(np.isnan(fitParams))):
			# theta = np.linspace(0, np.pi, 100)
			# plt.plot(theta * 180 / np.pi, VonMisses(theta * np.pi / 180, *fitParams), 'r')
			plt.title('neuron# %s, osi = %.4s, '%(i, osi) + r'$ \chi^2 = %.5s, \, SP = %s,$'%(chiSquare, int(IS_SINGLE_PEAK), ))
		    # plt.title('neuron# %s, osi = %.4s, '%(i, osi) + r'$ \chi^2 = %.5s$'%(chiSquare))
	    else:
		plt.title('neuron# %s, osi = %.4s, '%(i, osi))
            # ipdb.set_trace()		
	    plt.waitforbuttonpress()
	    plt.clf()
    if IF_FIT and not IF_PLOT:
	return tc, chiSquareArray, tcSem, tcSemOverAngles, avgAbsResidual
    else:
	return tc, tcSem, tcSemOverAngles

def ComputeTuningSEM2(p, gamma, nPhis, mExt, mExtOne, nChnks = 2, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_PLOT = False, IF_FIT = False, nNeurons = 10000, neuronType = 'E'):
    # plot tc with sem given atleast 2 chunks of simulations
    NE = N
    NI = N
    nNeurons = NE + NI
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    tc = np.empty((2, nNeurons, nPhis))
    tc[:] = np.nan
    for kChunk in range(nChnks):
	for i, iPhi in enumerate(phis):
	    fr = LoadFrChnk(p, gamma, iPhi, mExt, mExtOne, kChunk, trNo, T, NE, K, nPop)
	    tc[kChunk, :, i] = fr
    tcSem = np.squeeze(np.nanstd(tc, 0)) / np.sqrt(nChnks)
    tc = np.squeeze(np.nanmean(tc, 0))
    IS_GOOD_FIT = np.zeros((nNeurons, ))
    IS_SINGLE_PEAKED = np.zeros((nNeurons, ))
    if IF_FIT and not IF_PLOT:
        chiSquareArray = np.empty((nNeurons, ))
	chiSquareArray[:] = np.nan
	for k in range(nNeurons):
	    _, _, chiSquareArray[k], IS_GOOD_FIT[k] = FitVonMisses2(phis * np.pi / 180, tc[k, :], tcSem[k, :])
	    tmpTc = tc[k, :]
	    meanRate = np.nanmean(tc[k, :])
	    tmpTc = np.concatenate((tmpTc, [tmpTc[0], tmpTc[1]]))
	    peakIdx = argrelextrema(tmpTc, np.greater)[0]
	    if(peakIdx.size >= 2):
		secondPeak = np.sort(tmpTc[peakIdx])[-2]
		IS_SINGLE_PEAKED[k] = IS_GOOD_FIT[k] and ((np.nanmax(tc[k, :]) - meanRate) >= 0.05) and ((secondPeak < 0.05))
    if IF_PLOT:
        plt.ion()
    	for i in np.random.randint(0, NE, 101):
	    print 'neuron#', i
	    (_, caps, _) = plt.errorbar(np.concatenate((phis, [180])), np.concatenate((tc[i, :], [tc[i, 0]])), fmt = 'ko-', markersize = 3, yerr = np.concatenate((tcSem[i, :], [tcSem[i, 0]])), lw = 0.8, elinewidth=0.8)
	    for cap in caps:
		# cap.set_color('red')
		cap.set_markeredgewidth(0.8)
	    osi = OSI(tc[i, :], phis)
	    po = GetPhase(tc[i, :], phis, IF_IN_RANGE = True)
	    _, ymax = plt.ylim()
	    plt.vlines(po, 0, ymax, color = 'k')
            IS_SINGLE_PEAK = False
	    # ipdb.set_trace()
	    if IF_FIT:
		fitParams, _, chiSquare, IS_SINGLE_PEAK = FitVonMisses2(phis * np.pi / 180, tc[i, :], tcSem[i, :], IF_PLOT)
		if(~np.any(np.isnan(fitParams))):
		    meanRate = np.nanmean(tc[i, :])
		    
		    print chiSquare, np.nanmax(tc[i, :]) - meanRate, GammaQ(4, chiSquare), GammaQ(5, chiSquare)
		    tmpTc = tc[i, :]
		    tmpTc = np.concatenate((tmpTc, [tmpTc[0], tmpTc[1]]))
		    peakIdx = argrelextrema(tmpTc, np.greater)[0]
		    print 'neuron idx = ', i, 'peak idices = ', peakIdx
		    if(peakIdx.size >= 2):
			secondPeak = np.sort(tmpTc[peakIdx])[-2]
			print '2nd peak = ', secondPeak
			IS_SINGLE_PEAK = IS_SINGLE_PEAK and ((np.nanmax(tc[i, :]) - meanRate) >= 0.05) and ((secondPeak < 0.05))
			IS_SINGLE_PEAK2 = IS_SINGLE_PEAK and ((np.nanmax(tc[i, :]) - meanRate) >= 0.05) and ((secondPeak < 0.05))
		    if(~np.any(np.isnan(fitParams))):
			plt.title('neuron# %s, osi = %.4s, '%(i, osi) + r'$ \chi^2 = %.5s, \, SP = %s, sp2 = %s$'%(chiSquare, int(IS_SINGLE_PEAK), int(IS_SINGLE_PEAK2)))
	    else:
		plt.title('neuron# %s, osi = %.4s, '%(i, osi))
	    plt.waitforbuttonpress()
	    plt.clf()
    if IF_FIT and not IF_PLOT:
	return tc, tcSem, chiSquareArray, IS_SINGLE_PEAKED, IS_GOOD_FIT
    else:
	return tc, tcSem

def IF_SINGLE_PEAK(tcSemOverAngles, avgAbsResidual):
    peakyMetric = avgAbsResidual / tcSemOverAngles
    return peakyMetric < 1

def IF_SINGLE_PEAK2(chiSquare, tcSemOverAngles):
    peakyMetric = chiSquare / tcSemOverAngles
    return peakyMetric < 1

def PCA(X, ndims = 2, IF_PLOT = False):
    # X  : trial-by-dim
    centeredX = X - X.mean(0)
    covMat = np.dot(X.T, X)
    eigVal, eigVec = np.linalg.eig(covMat)
    projection = np.dot(X, eigVec[:, :ndims])
    if IF_PLOT:
        fig = plt.figure()
        if ndims == 2:
            plt.plot(projection[:, 0], projection[:, 1], 'k.')
        if ndims == 3:
            ax = fig.add_subplot(111, projection='3d')
            plt.plot(projection[:, 0], projection[:, 1], projection[:, 2], 'k.')
    return eigVal, eigVec, projection

# def 


def MultiPeaksDistr2(p, gamma, nPhis, mExt, mExtOne, nChnks = 2, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_PLOT = False, IF_FIT = True, peakRate = 0.10, IF_COMPUTE = False):
    if IF_COMPUTE:
	NE = N
	NI = N
	print 'computing tc...',
	sys.stdout.flush()
	tc, tcSem, chiSquared, IS_SINGLE_PEAK, IS_GOOD_FIT = ComputeTuningSEM2(p, gamma, nPhis, mExt, mExtOne, nChnks, trNo, N, K, nPop, T, IF_PLOT, IF_FIT)
	print ' done'
	responsiveNeruonIdx = np.nanmax(tc, 1) >= peakRate
        # IS_SINGLE_PEAK = np.logical_and(IS_SINGLE_PEAK, ((np.nanmax(tc, 1) - meanRate) >= 0.05))
	percentSinglePeakE = IS_SINGLE_PEAK[:NE].sum() * 100 / float(NE) #float(responsiveNeruonIdx[:NE].sum())
	percentSinglePeakI = IS_SINGLE_PEAK[NE:].sum() * 100 / float(NI) #float(responsiveNeruonIdx[NE:].sum())
	np.save('./data/p%sg%s_m0%s_mOne%s'%(int(10 * p), int(10 * gamma), int(1e3 * mExt), int(1e3 * mExtOne)), np.array([tc, chiSquared, tcSem, IS_SINGLE_PEAK, IS_GOOD_FIT, percentSinglePeakE, percentSinglePeakI]))
    else:
	tc, chiSquared, tcSem, IS_SINGLE_PEAK, IS_GOOD_FIT, percentSinglePeakE, percentSinglePeakI = np.load('./data/p%sg%s_m0%s_mOne%s.npy'%(int(10 * p), int(10 * gamma), int(1e3 * mExt), int(1e3 * mExtOne)))
    return percentSinglePeakE, percentSinglePeakI

def MultiPeaksDistr(p, gamma, nPhis, mExt, mExtOne, nChnks = 2, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_PLOT = False, IF_FIT = True, peakRate = 0.10, IF_COMPUTE = False):
    if IF_COMPUTE:
	NE = N
	NI = N
	print 'computing tc...',
	sys.stdout.flush()
	tc, chiSquared, tcSem, tcSemOverAngles, avgAbsResidual = ComputeTuningSEM(p, gamma, nPhis, mExt, mExtOne, nChnks, trNo, N, K, nPop, T, IF_PLOT, IF_FIT)
	print ' done'
	responsiveNeruonIdx = np.nanmax(tc, 1) >= peakRate
	nResponsiveNeurons = responsiveNeruonIdx.sum()


        print tcSemOverAngles.shape
	print avgAbsResidual

	nSinglePeakNeurons = 0
	# ipdb.set_trace()
	avgAbsRes = np.empty((NE + NI, ))
	for k in range(NE + NI):
	  avgAbsRes[k] = np.mean(avgAbsResidual[k] / tcSemOverAngles[k])
        meanRate = np.nanmean(tc, 1)
	
        IS_SINGLE_PEAK = np.logical_and(np.logical_or((avgAbsRes <= 2.0), (chiSquared <= 0.02)), ((np.nanmax(tc, 1) - meanRate) >= 0.05))


	
	# peakyMetric = np.sqrt(chiSquared) / tcSemOverAngles
	nSinglePeakNeuronsE = IS_SINGLE_PEAK.sum()
        #np.logical_and(responsiveNeruonIdx, peakyMetric < 1)
	percentSinglePeakE = IS_SINGLE_PEAK[:NE].sum() * 100 / float(responsiveNeruonIdx[:NE].sum())
	percentSinglePeakI = IS_SINGLE_PEAK[NE:].sum() * 100 / float(responsiveNeruonIdx[NE:].sum())
	np.save('./data/p%sg%s_m0%s_mOne%s'%(int(10 * p), int(10 * gamma), int(1e3 * mExt), int(1e3 * mExtOne)), np.array([tc, chiSquared, tcSem, tcSemOverAngles, avgAbsResidual, nSinglePeakNeurons, percentSinglePeakE, percentSinglePeakI]))
    else:
	tc, chiSquared, tcSem, tcSemOverAngles, avgAbsResidual, nSinglePeakNeurons, peakyMetric, percentSinglePeakE, percentSinglePeakI = np.load('./data/p%sg%s_m0%s_mOne%s.npy'%(int(10 * p), int(10 * gamma), int(1e3 * mExt), int(1e3 * mExtOne)))
    return percentSinglePeakE, nSinglePeakNeurons, peakyMetric

    
def ComputeTuning(p, gamma, nPhis, mExt, mExtOne, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_PLOT = False, IF_FIT = False):
    NE = N
    NI = N
    tc = np.zeros((NE + NI, nPhis))
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    for i, iPhi in enumerate(phis):
	fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, trNo, T, NE, K, nPop)
	# print fr.shape
	tc[:, i] = fr
    plt.ion()
    if IF_FIT and not IF_PLOT:
	chiSquareArray = np.empty((N, ))
	chiSquareArray[:] = np.nan
	for k in range(N):
	    _, _, chiSquareArray[k], _ = FitVonMisses(phis * np.pi / 180, tc[k, :])
    if IF_PLOT:
	for i in np.random.randint(0, NE, 101):
	    plt.plot(np.concatenate((phis, [180])), np.concatenate((tc[i, :], [tc[i, 0]])), 'ks-')
	    osi = OSI(tc[i, :], phis)
	    po = GetPhase(tc[i, :], phis, IF_IN_RANGE = True)
	    _, ymax = plt.ylim()
	    plt.vlines(po, 0, ymax, color = 'k')
	    if IF_FIT:
                fitParams, _, chiSquare, _ = FitVonMisses(phis * np.pi / 180, tc[i, :], IF_PLOT_FIT = True)
		if(~np.any(np.isnan(fitParams))):
		    # theta = np.linspace(0, np.pi, 100)
		    # plt.plot(theta * 180 / np.pi, VonMisses(theta * np.pi / 180, *fitParams), 'r')
		    plt.title('neuron# %s, osi = %.4s, '%(i, osi) + r'$ \chi^2 = %.5s$'%(chiSquare))
	    else:
		plt.title('neuron# %s, osi = %.4s, '%(i, osi))
	    plt.waitforbuttonpress()
	    plt.clf()
    if IF_FIT and not IF_PLOT:
	return tc, chiSquareArray
    else:
	return tc

def VonMisses(phi, baseLine, amplitude, width, po):
    return baseLine + amplitude * np.exp(width * (np.cos(2.0 * (phi - po))))

def FitVonMisses(x, y, ySEM, IF_PLOT_FIT = False):
    # 
    x = np.concatenate((x, [np.pi]))
    y = np.concatenate((y, [y[0]]))
    ySEM = np.concatenate((ySEM, [ySEM[0]]))
    
    # nParams = 4
    # degsOfFreedom = nPhis - nParams
    # significanceVal = 0.05
    # criticalChi2 = ChiSquared.isf(significanceVal, degsOfFreedom)
    try:
	bnds = ((0, 0, -np.inf, 0), (1, 1, np.inf, np.pi))
	#ipdb.set_trace()
	fitParams, fitError = curve_fit(VonMisses, x, y, bounds = bnds, max_nfev = 4000) #, maxfev = 4000)
	theta = np.linspace(0, np.pi, 100)
	fitY = VonMisses(x, *fitParams)
	chiSquare = ((y - fitY)**2 / ySEM).sum()
	# avgAbsResidual = np.mean(np.abs(y - fitY))
	avgAbsResidual = np.abs(y - fitY)
	# chiSquare = ((y - fitY)**2).sum()   
	if IF_PLOT_FIT:
	    plt.plot(theta * 180 / np.pi, VonMisses(theta, *fitParams), 'r')
    except RuntimeError, e:
	fitParams = np.empty((4, ))
	fitParams[:] = np.nan
	fitError = np.nan
	chiSquare = np.nan
	avgAbsResidual = np.empty((len(x) - 1, ))
	avgAbsResidual[:] = np.nan
	print e
    if(~np.any(np.isnan(fitParams))):
	return fitParams, fitError, chiSquare, avgAbsResidual[:-1]
    else:
	return fitParams, fitError, chiSquare, avgAbsResidual


def FitVonMisses2(x, y, ySEM, IF_PLOT_FIT = False):
    #
    x = np.concatenate((x, [np.pi]))
    y = np.concatenate((y, [y[0]]))
    ySEM = np.concatenate((ySEM, [ySEM[0]]))
    nPhis = len(x)
    nParams = 4
    degsOfFreedom = nPhis - nParams
    significanceVal = 0.05
    try:
	bnds = ((0, 0, -np.inf, 0), (1, 1, np.inf, np.pi))
	fitParams, fitError = curve_fit(VonMisses, x, y, bounds = bnds, max_nfev = 4000) #, maxfev = 4000)
	theta = np.linspace(0, np.pi, 100)
	fitY = VonMisses(x, *fitParams)
	vidx = y > 0
	chiSquare = ((y[vidx] - fitY[vidx])**2 / ySEM[vidx]).sum()
        qVal = GammaQ(degsOfFreedom, chiSquare) # probability the chi^2 is exceedes the computed value just by chance
	IS_GOOD_FIT = qVal > significanceVal 
	if IF_PLOT_FIT:
	    plt.ion()
	    plt.plot(theta * 180 / np.pi, VonMisses(theta, *fitParams), 'r')
	    plt.show()
    except RuntimeError, e:
	fitParams = np.empty((4, ))
	fitParams[:] = np.nan
	fitError = np.nan
	chiSquare = np.nan
	print e
    if(~np.any(np.isnan(fitParams))):
	return fitParams, fitError, chiSquare, IS_GOOD_FIT
    else:
	return fitParams, fitError, chiSquare, False

def ME1Analytic(JE0, JEE, JEI, gamma, p, m0One, mE0, mI0):
    num = -HPrime(mE0) * gamma * JEE * m0One
    alpha = JEE**2 * mE0 + JEI**2 * mI0
    denom = np.sqrt(alpha) * (1.0 - HPrime(mE0) * JEE * p / np.sqrt(alpha))
    return num / denom

def AlignTuningCurves(tc, nPhis = 8, NE = 10000, NI = 10000):
    prefferedOri = np.argmax(tc, 1)    
    tcMat = np.empty((NE + NI, nPhis))
    for kNeuron in np.arange(NE + NI):
	tcMat[kNeuron, :] = np.roll(tc[kNeuron, :], -1 * prefferedOri[kNeuron])
    return tcMat

def CenterTuningCurve(tc, nPhis = 8):
    maxFr = np.argmax(tc)
    tc = np.roll(tc, -1 * maxFr)
    return np.roll(tc, nPhis / 2)

def TieZeroPI(x, IF_DEGREES = False):
    theta = np.linspace(0, np.pi, len(x), endpoint = False)
    theta = np.concatenate((theta, [np.pi]))
    x = np.concatenate((x, [x[0]]))
    if IF_DEGREES:
        theta *= 180.0 / np.pi
    return theta, x

def PopAvgTuning(p, gamma, nPhis, mExt, mExtOne, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_NEW_FIG = True, color = 'k'):
    NE = N
    NI = N
    tc = np.zeros((NE + NI, nPhis))
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    for i, iPhi in enumerate(phis):
	fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, trNo, T, NE, K, nPop)
	print fr.shape
	tc[:, i] = fr
    prefferedOri = np.argmax(tc, 1)    
    tcMat = np.empty((NE + NI, nPhis))
    for kNeuron in np.arange(NE + NI):
	tcMat[kNeuron, :] = np.roll(tc[kNeuron, :], -1 * prefferedOri[kNeuron])
    meanE = np.nanmean(tcMat[:NE, :], 0)
    meanI = np.nanmean(tcMat[NE:, :], 0)
    rotateMeanBy = 4
    meanE = np.roll(meanE, rotateMeanBy)
    meanI = np.roll(meanI, rotateMeanBy)
    theta = np.linspace(-90, 90, nPhis, endpoint = False)
    print "mE = ", np.nanmean(meanE), "MI = ", np.nanmean(meanI)
    if IF_NEW_FIG:
	plt.figure()
    plt.plot(theta, meanE, 'o-', color = color, label = r"$m_0^{(1)} = %s, \, p = %s, \, \gamma = %s$"%(mExtOne, p, gamma))
    # plt.plot(theta, meanI, 'ro-', label = 'I')
    plt.title(r"$m_0 = %s$"%(mExt), fontsize = 16)
    plt.xlabel(r'$\phi$' + '(deg)')
    plt.ylabel(r'$m_E(\phi)$')
    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 10})    
    
def M1Component(x):
    out = np.nan
    if len(x) > 0:
	dPhi = np.pi / len(x)
	out = 2.0 * np.absolute(np.dot(x, np.exp(-2.0j * np.arange(len(x)) * dPhi))) / len(x)
  # zk = 2.0 * np.dot(firingRate[validIdx], np.exp(2j * atTheta[validIdx] * np.pi / 180)) / validIdx.sum()	
    # print 'm1 out: ', out
    return out

def SmoothedMeanM1(p, gamma, phis, mExt, mExtOne, NE = 10000, NI = 10000, K = 1000, nTr = 1, nPop = 2, T = 1000, neuronType = 'E', IF_M1_SQRD = False, IF_SMOOTHED = False):
    m1 = []
    winSize = int(float(NE) / 10.0)
    N = NE
    if neuronType == 'I':
        winSize = int(float(NI) / 10.0)
        N = NI
    nValidTrials = np.zeros((len(phis), ))
    print "--" * 25
    print "N = %s, K = %s, p = %s, m0 = %s, m0One = %s" %(N, K, p, mExt, mExtOne)
    print "trial#: ", 
    sys.stdout.flush()
    # phis = np.arange(0, 180.0, 22.5)
    # ipdb.set_trace()
    for n in range(nTr):
        print n
        m1OfPhi = []
        print 'PHI_0: ', 
        for phiIdx, atPhi in enumerate(phis):
            try:
                # print "LoadFr() ", p, gamma, atPhi, mExt, mExtOne, n, T, NE, K, nPop
                # print '#populations:', nPop
                # raise SystemExit
                m = LoadFr(p, gamma, atPhi, mExt, mExtOne, n, T, NE, K, nPop)
                # ipdb.set_trace()
                # print m
                if(neuronType == 'E'):
                    m = m[:N]
                else:
                    m = m[N:]
                if IF_SMOOTHED:
                    smoothedM = CirConvolve(m, winSize)
                    m1OfPhi.append(M1Component(smoothedM))
                else:
                    m1OfPhi.append(M1Component(m))                    
                # print m1OfPhi
                nValidTrials[phiIdx] += 1
                print atPhi, 
            except IOError:
                print 'x', #file not found'
            sys.stdout.flush()
        m1.append(np.nanmean(m1OfPhi))
    print "\n" + "--" * 25
    print 'mean m1 = ', np.mean(m1), 'SEM = ', np.nanstd(m1) / np.sqrt(float(nValidTrials[0]))
    if IF_M1_SQRD:
	print "sqrd"
	return np.nanmean(np.array(m1)**2), np.std(np.array(m1)**2), np.std(np.array(m1)**2) / np.sqrt(float(nValidTrials[0]))
    else:
	return np.nanmean(m1), np.nanstd(m1), np.nanstd(m1) / np.sqrt(float(nValidTrials[0]))
    
def M1vsp_Smoothed(pList, gamma, nPhis, mExt, mExtOne, NList = [10000], KList =[1000], trNo = 0, nPop = 2, T = 1000, nTrials = 1, IF_NONSMOOTHED = False, IF_NEW_FIG = True):
    m1 = []
    # plt.ion()
    if IF_NEW_FIG:
	plt.figure()
    IF_LEGEND = False
    nTr = nTrials
    # ipdb.set_trace()
    phis = np.linspace(0, 180.0, nPhis, endpoint = False)
    for k, kK in enumerate(KList):
        for n, nNE in enumerate(NList):
            m1 = []
            semM1 = []
            m1Sqrd = []
            semM1Sqrd = []
            for iip, p in enumerate(pList):
		if( not IF_NONSMOOTHED):
		    tmpM1, dummy, tmpsem = SmoothedMeanM1(p, gamma, phis, mExt, mExtOne, nNE, nNE, kK, nTr, nPop, T)
		    # tmpM1Sqrd, dummySqrd, tmpsemSqrd = SmoothedMeanM1(nNE, kK, p, mExt, nTr, T, IF_M1_SQRD = True)
		else:
		    tmpM1, dummy, tmpsem = NonSmoothedMeanM1(nNE, kK, p, mExt, nTr, T)
                m1.append(tmpM1)
                semM1.append(tmpsem)
                # print m1
                # m1Sqrd.append(tmpM1Sqrd)
                # semM1Sqrd.append(tmpsemSqrd)
                print '==' * 28; print '||' * 28; print '==' * 28
            m1 = np.array(m1)
	    # m1Sqrd = np.array(m1Sqrd)
            validPIdx = ~np.isnan(m1)
            semM1 = np.array(semM1)
            # semM1Sqrd = np.array(semM1Sqrd)
            # ipdb.set_trace()
            print m1[validPIdx].shape, np.sum(validPIdx)
            print pList[validPIdx]
            print semM1[validPIdx]
            if(np.sum(np.array(validPIdx)) > 0):
                IF_LEGEND = True
                (_, caps, _) = plt.errorbar(pList[validPIdx], m1[validPIdx], fmt = 'o-', markersize = 3, yerr = semM1[validPIdx], lw = 0.8, elinewidth=0.8, label = r'$N = %s, K = %s, \gamma = %s$'%(nNE, kK, gamma))
                for cap in caps:
                    # cap.set_color('red')
                    cap.set_markeredgewidth(0.8)
                outArray = np.empty((sum(validPIdx), 3))
                outArray[:] = np.nan
                outArray[:, 0] = pList[validPIdx]
                outArray[:, 1] = m1[validPIdx]
                outArray[:, 2] = semM1[validPIdx]                
                print outArray
                analysisDataFldr = '/homecentral/srao/Documents/code/binary/c/analysis/data/'
                if nPop == 1:
                    analysisDataFldr = analysisDataFldr + 'onepop/'
                if nPop == 2:
                    analysisDataFldr = analysisDataFldr + 'twopop/'
                np.save(analysisDataFldr +  "sim_mE1_vs_p_m0%s_N%s_K%s"%(int(mExt * 1e3), nNE, kK), outArray)
                # outArray2 = np.empty((sum(validPIdx), 3))                
		# outArray2[:, 0] = outArray[:, 0]
		# outArray2[:, 1] = m1Sqrd[validPIdx]
		# outArray2[:, 2] = semM1Sqrd[validPIdx]                
                # np.save("./data/analysisData/sim_mE1_sqrd_vs_p_m0%s_N%s_K%s"%(int(mExt * 1e3), nNE, kK), outArray2)		
            

    plt.xlabel(r'$p$', fontsize = 16)    
    plt.ylabel(r'$\left[ m^{(1)} \right]_{\phi}$', fontsize = 16) # + ' avg over 10 realizations')
    if IF_LEGEND:
        plt.legend(loc = 0, frameon = False, numpoints = 1)
    plt.xlim(0, pList[-1] + .5)
    # plt.ylim(0, mExt + 0.1)
    if IF_NONSMOOTHED:
	plt.title(r"$m_0 = %s, m_0^{(1)} = %s$"%(mExt, mExtOne), fontsize = 16)
        figname = 'non_smoothed_m1_m0%s'%(int(mExt * 1e3))	
    else:
	# plt.title(r"$m_0 = %s, smoothened$"%(mExt), fontsize = 16)
	plt.title(r"$m_0 = %s,\, m_0^{(1)} = %s, \,$"%(mExt, mExtOne) + 'smoothened', fontsize = 16)
        figname = 'smoothed_m1_m0%s'%(int(mExt * 1e3))	
    print 'saving as', figname
    figFolder = '/homecentral/srao/Documents/code/binary/c/analysis/figs/'
    if nPop == 1:
        figFolder = figFolder + 'onepop/'
    if nPop == 2:
        figFolder = figFolder + 'twopop/'
    #-------------------------------------#
    plt.savefig(figFolder + figname + '.png')
    plt.show()
    if IF_LEGEND:
        plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 10})
    plt.xlim(0, pList[-1] + .5)
    #plt.ylim(0, mExt + 0.1)
    ylims = plt.ylim()
    plt.yticks([0, ylims[1]])
    paperSize = [2.5, 2.0]
    figFormat = 'pdf'
    axPosition = [0.28, 0.23, .65, 0.65]
    print figFolder, figname
    Print2Pdf(plt.gcf(),  figFolder + figname,  paperSize, figFormat=figFormat, labelFontsize = 10, tickFontsize=8, titleSize = 10.0, IF_ADJUST_POSITION = True, axPosition = axPosition)


def CirConvolve(signal, windowSize):
    ker = np.concatenate((np.ones((windowSize, )), np.zeros((signal.size - windowSize, ))))
    return np.real(np.fft.ifft( np.fft.fft(signal)*np.fft.fft(ker))) * (1.0 / float(windowSize))
    


def PrintTuningBook(tc, neuronType, nNeurons, fname, NE = 10000, NI = 10000, nPhis = 8):
    doc = Document(fname)
    doc.packages.append(Package('geometry', options=['left=2cm', 'right=2cm']))

    nFigsPerPage = 10
    nPages = int(np.ceil(nNeurons / float(nFigsPerPage)))
    plt.figure()
    plt.ioff()
    theta = np.linspace(0, 180, nPhis, endpoint = False)
    doc = Document(fname)
    doc.packages.append(Package('geometry', options=['left=2cm', 'right=2cm']))
    width = r'.3\linewidth'
    if neuronType == 'E':
	rndNeurons = np.random.randint(0, NE, nNeurons)
    else:
	rndNeurons = np.random.randint(NE, NE + NI, nNeurons)	
    for kk in range(nPages):
	rndNeuronsx = rndNeurons[kk * nFigsPerPage: (kk + 1) * nFigsPerPage]
	print kk * nFigsPerPage, (kk + 1) * nFigsPerPage
	with doc.create(Figure(position='htbp')) as plot:
	    for i in rndNeuronsx:
		plt.plot(theta, tc[i, :], 'ko-')
		mainAxis = plt.gca()
		mainAxis.set_title('neuron#%s'%(i), fontsize = 10)
		with doc.create(SubFigure(position='b', width=NoEscape(width))) as figure:
		    figure.add_plot(width=NoEscape(r'\linewidth'), dpi = 300) #*args, **kwargs)
		plt.clf()
    
    doc.generate_pdf(clean_tex=False)

