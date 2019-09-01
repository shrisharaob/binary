# def ComputeFFInOutPOCorrParallelAux(kappa, p, gamma, nPhis, mExt, mExtOne, rewireType, N, K, nPop, T, NFF, JE0, JI0, IF_PLOT, trNo):
#     out = np.nan
#     try:
# 	theta = np.linspace(0, 180, nPhis, endpoint = False)
# 	uFF = rw.ComputeFFInput(nPhis, p, gamma, kappa, mExt, mExtOne, trNo, N, K, nPop, NFF, JE0, JI0, 0.2, rewireType, T) #nPhis, p, gamma, kappa, mExt, mExtOne, trNo, nPop, NFF, JE0, JI0, cFF, rewireType, T)
# 	# tcIn = GetInputTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa, 'FF')
# 	# ipdb.set_trace()
#         # ipdb.set_trace()	
#         poFF = POofPopulation(uFF, theta, IF_IN_RANGE = True) * np.pi / 180
# 	tcOut = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
# 	poOut = POofPopulation(tcOut, theta, IF_IN_RANGE = True) * np.pi / 180
# 	if IF_PLOT:
# 	    plt.plot(poFF[:N], poOut[:N], '.')
# 	    plt.show()
# 	    ipdb.set_trace()
# 	print 'computing ccc... ', trNo
# 	sys.stdout.flush()
# 	out = CircularCorrCoeff(poFF[:N], poOut[:N])
# 	print 'done', trNo, 'CCC=', out
#     except IOError:
# 	print ''
# 	return out
#     return out


import numpy as np
import pylab as plt



def ScatterPO(mExtOne=0.0375):
    tc, ff = LoadFFInput(mExtOne=mExtOne)
    poE = POofPopulation(tc[:10000])
    poFF = POofPopulation(ff[:10000])

    plt.plot(poFF, poE, 'k.')
    plt.savefig('./figs/PNAS/po_mExtOne%s.png'%(mExtOne))

    plt.show()
    

def LoadFFInput(mExtOne, trNo=1, N=10000, nPhis=8):
    fldr =  '/homecentral/srao/Documents/code/binary/c/twopop/PRX_data/data/rewire/N10000K1000/m075/mExtOne%s/kappa0/p0gamma0/T1/tr%s/'%(int(1e3 * mExtOne), int(trNo))

    #twopop/PRX_data/data/rewire/N10000K1000/m075/mExtOne%s/kappa0/p0gamma0/T1/tr0/'%(int(1e3 * mExtOne))
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    tc = np.zeros((2*N, nPhis))
    tcFF = np.zeros((2*N, nPhis))    
    for i, phi in enumerate(phis):
        filename = 'meanrates_theta%.6f_tr%s.txt'%(phi, trNo)
        tc[:, i] = np.loadtxt(fldr + filename)
        filename = 'meanFFinput_theta%.6f_tr%s_last.txt'%(phi, trNo)
        tcFF[:, i] = np.loadtxt(fldr + filename)        
    return tc, tcFF

def POofPopulation(tc, theta = np.arange(0.0, 180.0, 22.5), IF_IN_RANGE = False):
    # return value in degrees
    nNeurons, nAngles = tc.shape
    theta = np.linspace(0, 180, nAngles, endpoint = False)
    po = np.zeros((nNeurons, ))
    for kNeuron in np.arange(nNeurons):
        po[kNeuron] = GetPhase(tc[kNeuron, :], theta, IF_IN_RANGE)
    return po 



def GetPhase(firingRate, atTheta, IF_IN_RANGE = False):
    
    out = np.nan
    # zk = np.dot(firingRate, np.exp(2j * atTheta * np.pi / 180))
    # out = np.angle(zk) * 180.0 / np.pi
    # if IF_IN_RANGE:
    #     if(out < 0):
    #         out += 360
    # return out * 0.5
    
    thetas = np.arange(0, 180, 22.5) * np.pi / 180.0
    
    x = np.dot(np.cos(2 * thetas), firingRate)
    y = np.dot(np.sin(2 * thetas), firingRate)
    out = np.arctan2(y, x) * 180 / np.pi

    if x < 0 and y > 0:
         out = out / 2
    if x > 0 and y < 0:
        out = (out + 360) / 2
    if x < 0 and y < 0:
        out = (out + 360) / 2
    return out
