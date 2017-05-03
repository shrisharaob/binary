import numpy as np
import pylab as plt
import sys
from matplotlib.pyplot import cm


def GetFldrName(p, realizIdx, nPop = 1, N = 10000, K = 1000, m_ext = 0.15, gamma = 0, T = 1000):
    if nPop == 1:
	fldr = '/homecentral/srao/Documents/code/binary/c/onepop/'
    else:
        fldr = '/homecentral/srao/Documents/code/binary/c/twopop/'	
    fldr = fldr + "data/N%sK%sm0%sp%sgamma%sT%s/real%s/"%(N, K, int(m_ext * 1e3), int(p*10), int(gamma), int(T * 1e-3), int(realizIdx))
    return fldr

def PltM1(m1, ax0, ax1, lblTxt, N, pColor):
    m1Amp = m1[:, 0]
    m1Phase = m1[:, 1]
    intervalLength =  float(N / 1000.)
    xAxis = np.arange(0, m1Amp.size, dtype = 'float') / float(1000)
    ax0.plot(xAxis, m1Amp, label = lblTxt, lw = 0.5, c = pColor)
    ax1.plot(xAxis, m1Phase * 180 / np.pi, label = lblTxt, lw = 0.5, c = pColor)
    ax0.set_xlabel('time (avg updates/spin)')
    ax1.set_xlabel('time (avg updates/spin)')    
    ax0.set_ylabel(r'$m_I^{(1)}$')
    ax1.set_ylabel(r'$\Phi(deg)$')    
    # plt.show()

def PltM1Hist(m1, ax0, ax1, nPoints2Discard, lblTxt, nBins = 25):
    m1Amp = m1[nPoints2Discard:, 0]
    m1Phase = m1[nPoints2Discard:, 1]    
    ax0.hist(m1Amp, nBins, normed = 1, histtype = 'step', label = lblTxt)
    ax1.hist(m1Phase * 180 / np.pi, nBins, normed = 1, histtype = 'step', label = lblTxt)
    ax0.set_xlabel(r'$m_I^{(1)}$')
    ax1.set_xlabel(r'$\Phi(deg)$')    
    ax0.set_ylabel('density')
    ax1.set_ylabel('density')    
    # plt.show()

def LoadM1File(p, trialNo, realizIdx, nPop = 1, N = 10000, K = 1000, m_ext = 0.15, gamma = 0, T = 1000):
    fldr = GetFldrName(p, realizIdx, nPop, N, K, m_ext, gamma, T)
    # fldr = '/homecentral/srao/Documents/code/binary/c/onepop/'
    # fldr = fldr + "data/N%sK%sm0%sp%sgamma%sT%s/"%(N, K, int(m_ext * 1e3), int(p*10), int(gamma), int(T * 1e-3))    
    m1 = np.loadtxt(fldr + 'MI1_inst_tr%s.txt'%(trialNo), delimiter = ';')
    # print m1.shap
    return m1

def PltM1OverRealization(p, colors, m0 = 0.15, nTr = 4, nRealizations = 10, T = 1000):
    fg0, ax0 = plt.subplots()
    fg1, ax1 = plt.subplots()
    for j in range(nRealizations):
	print 'realiz# = ', j,
        sys.stdout.flush()	
	for i in range(nTr):
	    if i == 0:
		print ' tr# = ', i,
	    else:
		print i,
            sys.stdout.flush()			
	    m1 = LoadM1File(p, i, j, m_ext = m0, T = T)
	    PltM1(m1, ax0, ax1, 'realization#%s'%(j), 10000, colors[i])
	print ' '
#        plt.show()
    fldr = '/homecentral/srao/Documents/code/binary/c/onepop/'	
    ax0.set_title(r"$m_0=%s,p=%s, %s trials$"%(m0, p, nTr))
    ax1.set_title(r"$m_0=%s,p=%s, %s trials$"%(m0, p, nTr))
    # ax0.legend(loc = 0, frameon = False); plt.draw()
    # ax1.legend(loc = 0, frameon = False); plt.draw()	
    fg0.savefig(fldr + "figs/mi1_vs_time_p%s_m0%s_rzAll.png"%(int(p*10), int(m0*1e3)))
    fg1.savefig(fldr + "figs/mi1_phase_vs_time_p%s_m0%s_rzAll.png"%(int(p*10), int(m0*1e3)))
	# plt.figure(fg0.number); plt.clf()
	# plt.figure(fg1.number); plt.clf()

    plt.close('all')	


def PltM1OverTr(p, m0 = 0.15, nTr = 4, nRealizations = 10, T = 1000):
    for j in range(nRealizations):
	fg0, ax0 = plt.subplots()
	fg1, ax1 = plt.subplots()
	print 'realiz# = ', j,
        sys.stdout.flush()	
	for i in range(nTr):
	    if i == 0:
		print ' tr# = ', i,
	    else:
		print i,
            sys.stdout.flush()			
	    m1 = LoadM1File(p, i, j, m_ext = m0, T = T)
	    PltM1(m1, ax0, ax1, 'tr#%s'%(i))
	print ' '
	fldr = '/homecentral/srao/Documents/code/binary/c/onepop/'	
	ax0.set_title(r"$m_0=%s,p=%s,realization# = %s$"%(m0, p, j))
	ax1.set_title(r"$m_0=%s,p=%s,realization# = %s$"%(m0, p, j))
	ax0.legend(loc = 0, frameon = False); plt.draw()
        ax1.legend(loc = 0, frameon = False); plt.draw()	
	fg0.savefig(fldr + "figs/mi1_vs_time_p%s_m0%s_rz%s.png"%(int(p*10), int(m0*1e3), j))
	fg1.savefig(fldr + "figs/mi1_phase_vs_time_p%s_m0%s_rz%s.png"%(int(p*10), int(m0*1e3), j))
	# plt.figure(fg0.number); plt.clf()
	# plt.figure(fg1.number); plt.clf()	
        plt.close('all')	

	

# def PltM1Scatter(p, m0 = 0.15, nTr = 10):
#     fg0, ax0 = plt.subplots()
#     fg1, ax1 = plt.subplots()
#     m1 = LoadM1File(p, 0, m_ext = m0)    
#     for i in range(1, nTr):
# 	m2 = LoadM1File(p, i, m_ext = m0)
# 	ax0.scatter(m1[20000:, 0], m1[20000:, 0], '.')
# 	ax1.scatter(m1[20000:, 1], m1[20000:, 1], '.')	
#     fldr = '/homecentral/srao/Documents/code/binary/c/onepop/'	
#     ax0.set_title(r"$m_0 = %s, p = %s$"%(m0, p))
#     ax1.set_title(r"$m_0 = %s, p = %s$"%(m0, p))
#     fg0.savefig(fldr + "figs/mi1_scatter_p%s_m0%s.png"%(int(p*10), int(m0*1e3)))
#     fg1.savefig(fldr + "figs/mi1_phase_scatter_p%s_m0%s.png"%(int(p*10), int(m0*1e3)))    

def PltM1HistOverTr(p, m0 = 0.15, nTr = 4, nRealizations = 10, T = 1000):
    for j in range(nRealizations):
	fg0, ax0 = plt.subplots()
	fg1, ax1 = plt.subplots()
	print 'realiz# = ', j,
        sys.stdout.flush()	
	for i in range(nTr):
	    if i == 0:
		print ' tr# = ', i,
	    else:
		print i,
	    sys.stdout.flush()
	    m1 = LoadM1File(p, i, j, m_ext = m0, T = T)
	    PltM1Hist(m1, ax0, ax1, 10000, 'tr#%s'%(i))
        print ''     
	fldr = '/homecentral/srao/Documents/code/binary/c/onepop/'	
	ax0.set_title(r"$m_0 = %s, p=%s, realization# = %s$"%(m0, p, j))
	ax1.set_title(r"$m_0 = %s, p=%s, realization# = %s$"%(m0, p, j))
	ax0.legend(loc = 0, frameon = False); plt.draw()
        ax1.legend(loc = 0, frameon = False); plt.draw()	
	fg0.savefig(fldr + "figs/mi1_hist_p%s_m0%s_rz%s.png"%(int(p*10), int(m0*1e3), j))
	fg1.savefig(fldr + "figs/mi1_phase_hist_p%s_m0%s_rz%s.png"%(int(p*10), int(m0*1e3), j))
	# plt.figure(fg0.number); plt.clf()
	# plt.figure(fg1.number); plt.clf()	
        plt.close('all')


	

if __name__ == "__main__":
    p = float(sys.argv[1])
    T = int(sys.argv[2])
    nTr = int(sys.argv[3])
    color=cm.rainbow(np.linspace(0,1, 10))
    PltM1OverRealization(p, color, T = T, nTr = nTr, m0 = 0.15, nRealizations = 1)
    
    # PltM1HistOverTr(p = 3.5, m0 = 0.15, nTr = 2, nRealizations = 10)
    # PltM1OverTr(p = 3.5, m0 = 0.15, nTr = 2, nRealizations = 10)
