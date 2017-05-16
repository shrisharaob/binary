import numpy as np
import pylab as plt
import os, sys
import ipdb

basefolder = "/homecentral/srao/Documents/code/mypybox"
sys.path.append('/homecentral/srao/Documents/code/mypybox/utils')
from Print2Pdf import Print2Pdf


def LoadFr(p, gamma, phi, mExt, mExtOne, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2):
    baseFldr = '/homecentral/srao/Documents/code/binary/c/'
    if nPop == 1:
	baseFldr = baseFldr + 'onepop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if nPop == 2:
	baseFldr = baseFldr + 'twopop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    # print baseFldr
    filename = 'meanrates_theta%.6f_tr0.txt'%(phi)
    return np.loadtxt(baseFldr + filename)
 
def LoadSpkTimes(p, gamma, phi, mExt, mExtOne, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2):
    baseFldr = '/homecentral/srao/Documents/code/binary/c/'
    if nPop == 1:
	baseFldr = baseFldr + 'onepop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if nPop == 2:
	baseFldr = baseFldr + 'twopop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    print baseFldr
    return np.loadtxt(baseFldr + 'spktimes.txt', delimiter = ';')

    
def OSI(firingRate, atTheta):
    out = np.nan
    zk = np.dot(firingRate, np.exp(2j * atTheta * np.pi / 180))
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

def PltOSIHist(p, gamma, nPhis, mExt, mExtOne, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000):
    NE = N
    NI = N
    tc = np.zeros((NE + NI, nPhis))
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    for i, iPhi in enumerate(phis):
        print 'loading from fldr: ', 
	fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, trNo, T, NE, K, nPop)
        if(np.isnan(fr)):
            print 'file not found!'
	# print fr.shape
	tc[:, i] = fr
    osi = OSIOfPop(tc[:NE, :], phis)
    plt.figure()
    plt.hist(osi[~np.isnan(osi)], 27)
    plt.show()
    
def ComputeTuning(p, gamma, nPhis, mExt, mExtOne, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000):
    NE = N
    NI = N
    tc = np.zeros((NE + NI, nPhis))
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    for i, iPhi in enumerate(phis):
	fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, trNo, T, NE, K, nPop)
	# print fr.shape
	tc[:, i] = fr
    plt.ion()
    for i in np.random.randint(0, NE, 101):
	plt.plot(phis, tc[i, :], 'ks-')
	plt.waitforbuttonpress()
	plt.clf()

def M1Component(x):
    dPhi = np.pi / len(x)
    out = 2.0 * np.absolute(np.dot(x, np.exp(-2.0j * np.arange(len(x)) * dPhi))) / len(x)
    # print 'm1 out: ', out
    return out

def SmoothedMeanM1(p, gamma, nPhis, mExt, mExtOne, NE = 10000, NI = 10000, K = 1000, nTr = 1, nPop = 2, T = 1000, neuronType = 'E', IF_M1_SQRD = False, IF_SMOOTHED = False):
    m1 = []
    winSize = int(float(NE) / 10.0)
    N = NE
    if neuronType == 'I':
        winSize = int(float(NI) / 10.0)
        N = NI
    nValidTrials = np.zeros((nPhis, ))
    print "--" * 25
    print "N = %s, K = %s, p = %s, m0 = %s, m0One = %s" %(N, K, p, mExt, mExtOne)
    print "trial#: ", 
    sys.stdout.flush()
    phis = np.arange(0, 180.0, 22.5)
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
    
    print 'mean m1 = ', np.mean(m1), 'SEM = ', np.std(m1) / np.sqrt(float(nValidTrials[0]))
    if IF_M1_SQRD:
	print "sqrd"
	return np.nanmean(np.array(m1)**2), np.std(np.array(m1)**2), np.std(np.array(m1)**2) / np.sqrt(float(nValidTrials[0]))
    else:
	return np.nanmean(m1), np.std(m1), np.std(m1) / np.sqrt(float(nValidTrials[0]))

    
def M1vsp_Smoothed(pList, gamma, nPhis, mExt, mExtOne, NList, KList, trNo = 0, nPop = 2, T = 1000, nTrials = 1, IF_NONSMOOTHED = False):
    m1 = []
    # plt.ion()
    # plt.figure()
    IF_LEGEND = False
    nTr = nTrials
    # ipdb.set_trace()    
    for k, kK in enumerate(KList):
        for n, nNE in enumerate(NList):
            m1 = []
            semM1 = []
            m1Sqrd = []
            semM1Sqrd = []
            for iip, p in enumerate(pList):
		if( not IF_NONSMOOTHED):
		    tmpM1, dummy, tmpsem = SmoothedMeanM1(p, gamma, nPhis, mExt, mExtOne, nNE, nNE, kK, nTr, nPop, T)
		    # tmpM1Sqrd, dummySqrd, tmpsemSqrd = SmoothedMeanM1(nNE, kK, p, mExt, nTr, T, IF_M1_SQRD = True)
		else:
		    tmpM1, dummy, tmpsem = NonSmoothedMeanM1(nNE, kK, p, mExt, nTr, T)
                m1.append(tmpM1)
                semM1.append(tmpsem)
                print m1
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
                (_, caps, _) = plt.errorbar(pList[validPIdx], m1[validPIdx], fmt = 'o-', markersize = 3, yerr = semM1[validPIdx], lw = 0.8, elinewidth=0.8, label = 'N = %s, K = %s'%(nNE, kK))
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
    plt.ylabel(r'$\left[ m^{(1)} \right]_J$', fontsize = 16) # + ' avg over 10 realizations')
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
        plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 4})
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
    
