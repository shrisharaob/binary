import numpy as np
import pylab as plt
import sys
import Scripts as sc

IF_COMPUTE = int(sys.argv[1])
mExt = float(sys.argv[2])


nPhis = 8
NE = 10000

if IF_COMPUTE:
    mExtOne = float(sys.argv[3])
    p = float(sys.argv[4])
    gamma = float(sys.argv[5])
    out = sc.MultiPeaksDistr2(p, gamma, nPhis, mExt, mExtOne, IF_PLOT = 0, IF_FIT = 1, IF_COMPUTE = 1)
else:
    pList = range(9)
    gList = [0] #, 1, 2, 3, 4, 5, 6, 7] #range(7) #, 2, 4, 6, 7]
    mExtOneList = [0.0375, 0.75]
    fg0 = plt.figure()
    fg1 = plt.figure()
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

                        # responsiveIdx = np.max(tc, 1) > 0.1
			# nResponsiveE = responsiveIdx[:NE].sum()
			# nResponsiveI = responsiveIdx[NE:].sum()

			IS_SINGLE_PEAK_LOCAL = sc.GammaQ(5, chiSquared) > .05
			pcE = 100.0 * IS_SINGLE_PEAK_LOCAL[:NE].sum() / float(NE)
                	pcI = 100.0 * IS_SINGLE_PEAK_LOCAL[NE:].sum() / float(NE)


                        # IS_SINGLE_PEAK_LOCAL = np.logical_and(sc.GammaQ(4, chiSquared) > .05, responsiveIdx)
			# pcE = 100.0 * IS_SINGLE_PEAK_LOCAL[:NE].sum() / float(nResponsiveE)
                	# pcI = 100.0 * IS_SINGLE_PEAK_LOCAL[NE:].sum() / float(nResponsiveI)
			
			
			pcSinglePeakE.append(pcE)
			pcSinglePeakI.append(pcI)
			
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
                    plt.figure(fg0.number)
		    plt.plot(validPList, pcSinglePeakE, 'o-', label = r'$\gamma = %s$'%(gamma))
		    plt.figure(fg1.number)
		    plt.plot(validPList, pcSinglePeakI, 'o-', label = r'$\gamma = %s$'%(gamma))
		    # plt.plot(validPList, pcSinglePeakI, 's-', label = r'$I, \, m_0^{(1)} = %s, \gamma = %s$'%(mExtOne, gamma))

	    plt.figure(fg0.number)
	    plt.title(r"$m_0^{(0)} = %s,\, m_0^{(1)} = %s$"%(mExt, mExtOne), fontsize = 18)
	    plt.xlabel('p', fontsize = 18)
	    plt.ylabel('Single peaked(%), E', fontsize = 18)
	    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 12})
	    plt.gca().set_position([0.2, 0.2, .65, .65])    
	    plt.savefig('./figs/twopop/SinglePeak_Summary_E_m0%s_mOne%s.png'%(mExt, mExtOne))
	    
	    plt.figure(fg1.number)
   	    plt.title(r"$m_0^{(0)} = %s,\, m_0^{(1)} = %s$"%(mExt, mExtOne), fontsize = 18)
	    plt.xlabel('p', fontsize = 18)
	    plt.ylabel('Single peaked(%), I', fontsize = 18)
	    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 12})
	    plt.gca().set_position([0.2, 0.2, .65, .65])    
	    plt.savefig('./figs/twopop/SinglePeak_Summary_I_m0%s_mOne%s.png'%(mExt, mExtOne))
	    plt.show()
    except IOError, errorMsg:
	print "hello there: ", errorMsg
    
