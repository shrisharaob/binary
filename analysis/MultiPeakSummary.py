import numpy as np
import pylab as plt
import sys
import Scripts as sc
import FitGaussian as fg


def MultiPeaksDistr(p, gamma, nPhis, mExt, mExtOne, nChnks = 10, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_PLOT = False, IF_FIT = True, peakRate = 0.10, IF_COMPUTE = False, IF_COMPUTE_PC = False):

    NE = N
    NI = N
    if IF_COMPUTE:
	print 'computing tc...',
	sys.stdout.flush()
	# tc, tcSem, chiSquared, IS_SINGLE_PEAK, IS_GOOD_FIT = ComputeTuningSEM3(p, gamma, nPhis, mExt, mExtOne, nChnks, trNo, N, K, nPop, T, IF_PLOT, IF_FIT)
        # tc, tcSem, chiSquareArray, IS_GOOD_FIT, IS_RESPONSIVE
        tc, tcSem, chiSquareArray, IS_GOOD_FIT, IS_RESPONSIVE = fg.ComputeTuningSEM(p, gamma, nPhis, mExt, mExtOne, IF_PLOT = 0, IF_FIT = 1, trNo = trNo, T = T, nChnks = nChnks)
	print ' done'
	percentSinglePeakE = IS_RESPONSIVE[:NE].sum() * 100 / float(NE)
	percentSinglePeakI = IS_RESPONSIVE[NE:].sum() * 100 / float(NI)
	# np.save('./data/p%sg%s_m0%s_mOne%s'%(int(10 * p), int(10 * gamma), int(1e3 * mExt), int(1e3 * mExtOne)), np.array([tc, chiSquared, tcSem, IS_SINGLE_PEAK, IS_GOOD_FIT, percentSinglePeakE, percentSinglePeakI]))
    if IF_COMPUTE_PC:
	out = np.load('./data/fit_p%sg%s_m0%s_mOne%s_nPhis%s.npz'%(int(10 * p), int(10 * gamma), int(1e3 * mExt), int(1e3 * mExtOne), nPhis))
        tc = out['tc']
        tcSem = out['tcSem']
        IS_RESPONSIVE = out['IS_RESPONSIVE']
        IS_GOOD_FIT = out['IS_GOOD_FIT']
	fitParams = out['fitParams']
	thetas = np.linspace(0, 180, nPhis, endpoint = False)
	chiSquare = []
	fitError = []
	for i in range(NE):
	    fity = fg.PeriodicGaussian(thetas, *fitParams[i, :])
	    chiSquare.append(ChiSquare(tc[i, :], tcSem[i, :], fity))
	    fitMax = np.max(fg.PeriodicGaussian(np.linspace(0, 180, 100, endpoint = 0), *fitParams[i, :]))
	    fitError.append(np.sqrt(np.mean((tc[i, :] - fity) **2)) / fitMax)
	    tcMax = np.max(tc[i, :])
	    peakError = 100.0 * np.abs(tcMax - fitMax) / tcMax
	    plt.ion()
	    if fitError[i] < 0.06 and peakError < 5:
		plt.plot(thetas, tc[i, :], 'k.-')
		plt.plot(thetas, fity, 'go-')
		print fitError[i]
		# print chiSquare[i], fg.GammaQ(nPhis - 4, chiSquare[i])
		plt.waitforbuttonpress()
		plt.clf()
        percentSinglePeakE = float(np.sum((fitError[:NE] < 0.06))) / float(NE)
	percentSinglePeakI = float(np.sum((fitError[NE:] < 0.06))) / float(NI)
    else:
        out = np.load('./data/fit_p%sg%s_m0%s_mOne%s_nPhis%s.npz'%(int(10 * p), int(10 * gamma), int(1e3 * mExt), int(1e3 * mExtOne), nPhis))
        tc = out['tc']
        tcSem = out['tcSem']
        IS_RESPONSIVE = out['IS_RESPONSIVE']
        IS_GOOD_FIT = out['IS_GOOD_FIT']
	chiSquare = out['chiSquareArray']
	fitParams = out['fitParams']
	thetas = np.linspace(0, 180, nPhis, endpoint = False)
	chiSquare = []
	fitError = []
	peakError = []
	for i in range(NE + NI):
	    fity = fg.PeriodicGaussian(thetas, *fitParams[i, :])
	    fitMax = np.max(fg.PeriodicGaussian(np.linspace(0, 180, 100, endpoint = 0), *fitParams[i, :]))
	    tcMax = np.max(tc[i, :])	    
	    fitError.append(np.sqrt(np.mean((tc[i, :] - fity) **2)) / tcMax)

	    peakError.append(100.0 * np.abs(tcMax - fitMax) / tcMax)
	# fg.ipdb.set_trace()
	peakError = np.array(peakError)
        fitError = np.array(fitError)
	# errorMeasure = np.logical_and(fitError < 0.06, peakError < 5)
	errorMeasure = fitError < 0.1
	distantMeasure = fitError[errorMeasure[:NE]]
	print 'distantMeasure shape', distantMeasure.shape
        percentSinglePeakE = np.mean((errorMeasure[:NE]))
	percentSinglePeakI = np.mean((errorMeasure[NE:]))
	print percentSinglePeakE, percentSinglePeakI
    return percentSinglePeakE * 100, percentSinglePeakI * 100, fitError, distantMeasure

def PlotCumulativeDistr(x, IF_NEW_FIGURE = True, label = ''):
    x = np.sort(x)
    yvals = np.arange(1, len(x) + 1)/ float(len(x))
    if IF_NEW_FIGURE:
	plt.figure()
    plt.plot(x, yvals, label = label)
    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 10})
    plt.xlabel('dist')
    plt.ylabel('prob')

def ChiSquare(y, ysem, fity):
    vidx = ysem > 0
    ym = y[vidx]# - y[vidx].mean()
    fitym = fity[vidx]# - fity[vidx].mean()
    yMean = y[vidx].mean()
    fyMean = fity[vidx].mean()
    # return np.sqrt((ym - fitym) **2) / np.max(
    # return (((ym - fitym)**2/ ysem[vidx])).sum()
    return (((ym - fitym)**2/ fitym)).sum()
    # return (((y[vidx] - fity[vidx])**2 / ysem[vidx])).sum()
    # return ((y[vidx] - fity[vidx])**2).sum() # / ysem[vidx]).sum()

if __name__ == "__main__":
    IF_COMPUTE = int(sys.argv[1])
    mExt = float(sys.argv[2])
    nPhis = 36
    NE = 10000
    if IF_COMPUTE:
	mExtOne = float(sys.argv[3])
	p = float(sys.argv[4])
	gamma = float(sys.argv[5])
	trialNo = int(sys.argv[6])
	T = int(sys.argv[7])
	nChnks = int(sys.argv[8])
	out = MultiPeaksDistr(p, gamma, nPhis, mExt, mExtOne, IF_PLOT = 0, IF_FIT = 1, IF_COMPUTE = 1, trNo = trialNo, T = T, nChnks = nChnks)
    else:
	pList = [0, 7] #range(9)
	gList = [0] #, 1, 2, 3, 4, 5, 6, 7] #range(7) #, 2, 4, 6, 7]
	mExtOneList = [0.075] #, 0.075]
	fg0 = plt.figure()
	fg1 = plt.figure()
	IF_CDF_FIGURE = False
	try:
	    for mExtOne in mExtOneList:
		for gamma in gList:
		    pcSinglePeakE = []
		    pcSinglePeakI = []
		    validPList = []
		    validGList = []
		    for p in pList:
			print 'gamma= ', gamma, ' p=', p, ' <ChiSquare> = ', 
			try:
			    pcE, pcI, chiSquare, distantMeasure = MultiPeaksDistr(p, gamma, nPhis, mExt, mExtOne, IF_PLOT = 0, IF_FIT = 0, IF_COMPUTE = 0, IF_COMPUTE_PC = 0)
			    print np.nanmean(chiSquare)
			    # plt.figure()
			    # plt.hist(chiSquare, 100)
			    # plt.title(r"$m_0^{(0)} = %s,\, m_0^{(1)} = %s, p= %s, \gamma = %s$"%(mExt, mExtOne, p, gamma), fontsize = 18)
			    pcSinglePeakE.append(pcE)
			    pcSinglePeakI.append(pcI)
			    validPList.append(p)
			    validGList.append(gamma)

                            label = r'$p = %s, \gamma = %s$'%(p, gamma)
			    if not IF_CDF_FIGURE:
				PlotCumulativeDistr(distantMeasure, not IF_CDF_FIGURE, label = label)
				IF_CDF_FIGURE = True
			    else:
				PlotCumulativeDistr(distantMeasure, False, label = label)
                        
			except IOError, errorMsg:
			    print"hello there: ", errorMsg
		    print validPList
		    print 'pc E = ', pcSinglePeakE
		    print pcSinglePeakI	    
		    if(len(validPList) > 0):
			plt.figure(fg0.number)
			plt.plot(validPList, pcSinglePeakE, 'o-', label = r'$\gamma = %s$'%(gamma))
			plt.figure(fg1.number)
			plt.plot(validPList, pcSinglePeakI, 'o-', label = r'$\gamma = %s$'%(gamma))
			# plt.plot(validPList, pcSinglePeakI, 's-', label = r'$I, \, m_0^{(1)} = %s, \gamma = %s$'%(mExtOne, gamma))
		plt.figure(fg0.number)
		# plt.title(r"$m_0^{(0)} = %s,\, m_0^{(1)} = %s$"%(mExt, mExtOne), fontsize = 18)
		plt.xlabel('p', fontsize = 18)
		plt.ylabel('Single peaked(%), E', fontsize = 18)
		plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 12})
		plt.gca().set_position([0.2, 0.2, .65, .65])    
		plt.savefig('./figs/twopop/SinglePeak_Summary_E_m0%s_mOne%s.png'%(mExt, mExtOne))
		plt.figure(fg1.number)
		# plt.title(r"$m_0^{(0)} = %s,\, m_0^{(1)} = %s$"%(mExt, mExtOne), fontsize = 18)
		plt.xlabel('p', fontsize = 18)
		plt.ylabel('Single peaked(%), I', fontsize = 18)
		plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 12})
		plt.gca().set_position([0.2, 0.2, .65, .65])    
		plt.savefig('./figs/twopop/SinglePeak_Summary_I_m0%s_mOne%s.png'%(mExt, mExtOne))
		plt.show()
	except IOError, errorMsg:
	    print "hello there: ", errorMsg
    
