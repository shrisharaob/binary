import numpy as np
import pylab as plt
import sys
import Scripts as sc


mExt = float(sys.argv[1])
mExtOne = float(sys.argv[2])
p = float(sys.argv[3])
gamma = float(sys.argv[4])
IF_COMPUTE = int(sys.argv[5])
nPhis = 8
NE = 10000

if IF_COMPUTE:
    out = sc.MultiPeaksDistr2(p, gamma, nPhis, mExt, mExtOne, IF_PLOT = 0, IF_FIT = 1, IF_COMPUTE = 1)
else:
    pList = range(9)
    gList = [0, 1, 2, 3, 4, 5, 6, 7] #range(7) #, 2, 4, 6, 7]
    mExtOneList = [0.075]
    try:
	for mExtOne in mExtOneList:
	    for gamma in gList:
		pcSinglePeakE = []
		pcSinglePeakI = []
		validPList = []
		validGList = []
		for p in pList:
		    print 'gamma= ', gamma, ' p=', p
		    try:
			tc, chiSquared, tcSem, IS_SINGLE_PEAK, IS_GOOD_FIT, percentSinglePeakE, percentSinglePeakI = np.load('./data/p%sg%s_m0%s_mOne%s.npy'%(int(10 * p), int(10 * gamma), int(1e3 * mExt), int(1e3 * mExtOne)))

			IS_SINGLE_PEAK_LOCAL = sc.GammaQ(5, chiSquared) > .05
			pcE = 100.0 * IS_SINGLE_PEAK_LOCAL[:NE].sum() / float(NE)
                	pcI = 100.0 * IS_SINGLE_PEAK_LOCAL[NE:].sum() / float(NE)
			pcSinglePeakE.append(pcI)
			
			# pcSinglePeakE.append(percentSinglePeakE)
			# pcSinglePeakI.append(percentSinglePeakI)
			validPList.append(p)
			validGList.append(gamma)
		    except IOError, errorMsg:
			print"hello there: ", errorMsg
		print validPList
		print pcSinglePeakE
		print pcSinglePeakI	    
		if(len(validPList) > 0):
		    plt.plot(validPList, pcSinglePeakE, 'o-', label = r'$E, \, m_0^{(1)} = %s, \gamma = %s$'%(mExtOne, gamma))
		    # plt.plot(validPList, pcSinglePeakI, 's-', label = r'$I, \, m_0^{(1)} = %s, \gamma = %s$'%(mExtOne, gamma))

    except IOError, errorMsg:
	print "hello there: ", errorMsg


    plt.xlabel('p')
    plt.ylabel('Single peaked(%)')
    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 10})
    plt.savefig('./figs/SinglePeak_Summary.png')
    plt.show()	
    
