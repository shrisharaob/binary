import sys
from scipy.optimize import fsolve
import numpy as np
from scipy.integrate import quad, simps
from scipy.special import erfc, erfcinv
import pylab as plt
import ipdb
basefolder = "/homecentral/srao/Documents/code/mypybox"
sys.path.append(basefolder)
sys.path.append(basefolder + "/utils")
from Print2Pdf import Print2Pdf
import Scripts as sc
from multiprocessing import Pool

def MyErfc(z):
    # normal CDF
    return 0.5 * erfc(z / np.sqrt(2.0))

def MyErfc2(z):
    # normal CDF sqrd
    return (MyErfc(z)) **2    

def HPrime(z):
    return -1.0 * np.exp(-0.50 * z**2) / np.sqrt(2.0 * np.pi)

def QFuncInv(z):
    return np.sqrt(2.0) * erfcinv(z * 2.0)

def PrintError(outStruct, varName):
    print "-----------------------------------"
    print "Variable        :",  varName
    print "Solution        :", outStruct[0]
    print "error           :", outStruct[1]['fvec']
    print "Solution Status :", outStruct[-1]
    print "-----------------------------------"    

def CheckBalConditions(JEE, JEI, JIE, JII, JE0, JI0):
    JE = -JEI/JEE
    JI = -JII/JIE
    E = JE0
    I = JI0
    if((JE < JI) or (E/I < JE/JI) or (JE < 1) or (E/I < 1) or (JE/JE < 1)):
        print "NOT IN BALANCED REGIME!!!!!! "
        raise SystemExit

def alphaF(phi, EorI):
    out = np.dot(Jab**2, mFourier0 + p[EorI] * mFourier1 * np.cos(2 * phi))
    out = out[EorI]
    out += cFF * Ea[EorI]**2 * mExtZero     # add FF term
    return np.abs(out)

# def betaOfPhi(phi, q0, q1, EorI):
#     out = np.dot(Jab **2, q0 + p[EorI] * q1 * np.cos(2 * phi))
#     qExtZero = mExtZero**2 + 0.5 * mExtOne**2 # i.e qExt(Delta = 0)
#     qExtOne = 2 * mExtZero * mExtOne
#     out = out[EorI]    
#     out += cFF * Ea[EorI]**2 * (qExtZero + gamma[EorI] * qExtOne * np.cos(2 * phi))
#     return np.abs(out)

def betaOfPhi(phi, delta, q0, q1, EorI):
    out = np.dot(Jab **2, q0 + p * q1 * np.cos(2 * (phi + delta)))
    qExtZero = mExtZero**2 + 0.5 * mExtOne**2 * np.cos(4 * delta)
    qExtOne = 2 * mExtZero * mExtOne * np.cos(4 * delta)
    out = out[EorI]    
    out += cFF * Ea[EorI]**2 * (qExtZero + gamma[EorI] * qExtOne * np.cos(2 * phi))
    return np.abs(out)

def meanInput(uintial, *args):
    u0, u1 = uintial
    mFourier0, mFourier1 = args
    func0 = lambda phi: ((1.0 / (np.pi)) * MyErfc((1.0 - u0 - u1 * np.cos(2 * phi)) / np.sqrt(alphaF(phi, EorI = 0))))
    func1 = lambda phi: ((2.0 / (np.pi)) * MyErfc((1.0 - u0 - u1 * np.cos(2 * phi)) / np.sqrt(alphaF(phi, EorI = 1))) * np.cos(2 * phi))
    a = mFourier0 - quad(func0, 0, np.pi)[0]
    b = mFourier1 - quad(func1, 0, np.pi)[0]
    return (a, b)

def QuenchedDisorderOLD(qinitial, *args):
    q0E, q0I, q1E, q1I = qinitial
    mFourier0, mFourier1, u0E, u0I, u1E, u1I, delta = args
    q0 = np.array([q0E, q0I])
    q1 = np.array([q1E, q1I])
    hermiteDeg = 30
    hermiteX, hermiteWeights = np.polynomial.hermite.hermgauss(hermiteDeg)
    hermiteX = hermiteX * np.sqrt(2.0)
    hermiteWeights = hermiteWeights / np.sqrt(np.pi)
    
    uFuncE = lambda phi, delta: (u0E + u1E * np.cos(2 * (phi + delta)))
    uFuncI = lambda phi, delta: (u0I + u1I * np.cos(2 * (phi + delta)))    

    func0E = lambda phi: (1.0 / np.pi) * np.dot(MyErfc((1 - uFuncE(phi, delta) - np.sqrt(betaOfPhi(phi, delta, q0, q1, EorI = 0)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, delta, q0, q1, EorI = 0))) * MyErfc((1 - uFuncE(phi, -delta) - np.sqrt(betaOfPhi(phi, -delta, q0, q1, EorI = 0)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, -delta, q0, q1, EorI = 0))), hermiteWeights)
    
    func0I = lambda phi: (1.0 / np.pi) * np.dot(MyErfc((1 - uFuncI(phi, delta) - np.sqrt(betaOfPhi(phi, delta, q0, q1, EorI = 1)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, delta, q0, q1, EorI = 1))) * MyErfc((1 - uFuncI(phi, -delta) - np.sqrt(betaOfPhi(phi, -delta, q0, q1, EorI = 1)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, -delta, q0, q1, EorI = 1))), hermiteWeights)

    func1E = lambda phi: (2.0 / np.pi) * np.dot(MyErfc((1 - uFuncE(phi, delta) - np.sqrt(betaOfPhi(phi, delta, q0, q1, EorI = 0)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, delta, q0, q1, EorI = 0))) * MyErfc((1 - uFuncE(phi, -delta) - np.sqrt(betaOfPhi(phi, -delta, q0, q1, EorI = 0)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, -delta, q0, q1, EorI = 0))) * np.cos(2 * phi), hermiteWeights)
    
    func1I = lambda phi: (2.0 / np.pi) * np.dot(MyErfc((1 - uFuncI(phi, delta) - np.sqrt(betaOfPhi(phi, delta, q0, q1, EorI = 1)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, delta, q0, q1, EorI = 1))) * MyErfc((1 - uFuncI(phi, -delta) - np.sqrt(betaOfPhi(phi, -delta, q0, q1, EorI = 1)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, -delta, q0, q1, EorI = 1)))  * np.cos(2 * phi), hermiteWeights) #* np.cos(2 * phi)
    aE = quad(func0E, 0, np.pi)[0] - q0E
    aI = quad(func0I, 0, np.pi)[0] - q0I
    bE = quad(func1E, 0, np.pi)[0] - q1E
    bI = quad(func1I, 0, np.pi)[0] - q1I
    return (aE, aI, bE, bI)


def QuenchedDisorder(qinitial, *args):
    q0E, q0I, q1E, q1I = qinitial
    mFourier0, mFourier1, u0E, u0I, u1E, u1I, delta = args
    q0 = np.array([q0E, q0I])
    q1 = np.array([q1E, q1I])
    hermiteDeg = 30
    hermiteX, hermiteWeights = np.polynomial.hermite.hermgauss(hermiteDeg)
    hermiteX = hermiteX * np.sqrt(2.0)
    hermiteWeights = hermiteWeights / np.sqrt(np.pi)
    
    uFuncE = lambda phi, delta: (u0E + u1E * np.cos(2 * (phi + delta)))
    uFuncI = lambda phi, delta: (u0I + u1I * np.cos(2 * (phi + delta)))    

    func0E = lambda phi: (1.0 / np.pi) * np.dot(MyErfc((1 - uFuncE(phi, delta) - np.sqrt(betaOfPhi(phi, delta, q0, q1, EorI = 0)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, delta, q0, q1, EorI = 0))) * MyErfc((1 - uFuncE(phi, -delta) - np.sqrt(betaOfPhi(phi, delta, q0, q1, EorI = 0)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, delta, q0, q1, EorI = 0))), hermiteWeights)
    
    func0I = lambda phi: (1.0 / np.pi) * np.dot(MyErfc((1 - uFuncI(phi, delta) - np.sqrt(betaOfPhi(phi, delta, q0, q1, EorI = 1)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, delta, q0, q1, EorI = 1))) * MyErfc((1 - uFuncI(phi, -delta) - np.sqrt(betaOfPhi(phi, delta, q0, q1, EorI = 1)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, delta, q0, q1, EorI = 1))), hermiteWeights)

    func1E = lambda phi: (2.0 / np.pi) * np.dot(MyErfc((1 - uFuncE(phi, delta) - np.sqrt(betaOfPhi(phi, delta, q0, q1, EorI = 0)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, delta, q0, q1, EorI = 0))) * MyErfc((1 - uFuncE(phi, -delta) - np.sqrt(betaOfPhi(phi, delta, q0, q1, EorI = 0)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, delta, q0, q1, EorI = 0))) * np.cos(2 * phi), hermiteWeights)
    
    func1I = lambda phi: (2.0 / np.pi) * np.dot(MyErfc((1 - uFuncI(phi, delta) - np.sqrt(betaOfPhi(phi, delta, q0, q1, EorI = 1)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, delta, q0, q1, EorI = 1))) * MyErfc((1 - uFuncI(phi, -delta) - np.sqrt(betaOfPhi(phi, delta, q0, q1, EorI = 1)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, delta, q0, q1, EorI = 1)))  * np.cos(2 * phi), hermiteWeights) #* np.cos(2 * phi)
    aE = quad(func0E, 0, np.pi)[0] - q0E
    aI = quad(func0I, 0, np.pi)[0] - q0I
    bE = quad(func1E, 0, np.pi)[0] - q1E
    bI = quad(func1I, 0, np.pi)[0] - q1I
    return (aE, aI, bE, bI)

def FSolveAtDelta2(atDelta):
    # Returns q(0)OfDelta and q(1)OfDelta
    # print atDelta
    u0E, u1E =  fsolve(meanInput, uInitialGuess, args = (mFourier0[0], mFourier1[0]))
    u0I, u1I =  fsolve(meanInput, uInitialGuess, args = (mFourier0[1], mFourier1[1]))
    q0E, q0I, q1E, q1I =  fsolve(QuenchedDisorder, (0.01, 0.02, 0.00, 0.00), args = (mFourier1, mFourier1, u0E, u0I, u1E, u1I, atDelta))
    q0 = np.array([q0E, q0I])
    q1 = np.array([q1E, q1I])
    return  (q0, q1)

def GenTuningCurves(u0, u1, alpha, betaAtPhi, nNeurons, nPhis = 36):
    tc = np.zeros((nNeurons, nPhis))
    thetas = np.linspace(0, np.pi, nPhis, endpoint = False)
    xi = np.random.randn(nNeurons)
    plt.ion()
    plt.show()
    for i in range(nNeurons):
	for j, iPhi in enumerate(thetas):
	    u = u0 + u1 * np.cos(2 * iPhi)
	    beta = betaAtPhi(iPhi)[0]
	    tc[i, j] = MyErfc((1 - u - np.sqrt(beta) * xi[i]) / np.sqrt(alpha - beta))
	plt.plot(thetas, tc[i, :])
	plt.waitforbuttonpress()
	plt.clf()
    return tc

def TuningCurve(nNeurons): # = 1000, IF_PLOT = False):
    IF_PLOT = False
    u0E, u1E =  fsolve(meanInput, uInitialGuess, args = (mFourier0[0], mFourier1[0]))
    u0I, u1I =  fsolve(meanInput, uInitialGuess, args = (mFourier0[1], mFourier1[1]))
    q0E, q0I, q1E, q1I =  fsolve(QuenchedDisorder, (0.01, 0.02, 0.00, 0.00), args = (mFourier1, mFourier1, u0E, u0I, u1E, u1I, 0))
    q0 = np.array([q0E, q0I])
    q1 = np.array([q1E, q1I])    
    betaAtDeltaZero = lambda atPhi: betaOfPhi(atPhi, 0, q0, q1, EorI = 0)
    miOfPhi = np.zeros((atPhis.size, ))
    osi = []
    for l in range(nNeurons):
        qNoise = GenQuenchedNoiseWithCorr(maxCorrM)        
        for i, iPhi in enumerate(atPhis):
            miOfPhi[i] = MyErfc((1.0 - u0E - u1E * np.cos(2 * iPhi) - qNoise(iPhi)) / np.sqrt(alphaF(iPhi, 0) - betaAtDeltaZero(iPhi)))
	if IF_PLOT:
            plt.ion()	    
	    plt.plot(atPhis * 180 / np.pi, miOfPhi, 'ko-')
	    plt.grid()
	    plt.xlabel(r'$\phi$(deg)', fontsize = 20)
	    plt.ylabel(r'$m^i_E(\phi)$', fontsize = 20)
	    plt.waitforbuttonpress()
	    plt.clf()
	osi.append(OSI(miOfPhi, atPhis))
    return osi

def QMN(m, qnOfDelta, EorI):
    # qnOfDelta must be 2-by-nDelta matrix
    qn0 = (2. / np.pi) * simps(qnOfDelta[0, :] * np.cos(2 * m * atDeltas), atDeltas)
    qn1 = (p[EorI] * 2.0 / np.pi) * simps(qnOfDelta[1, :] * np.cos(2 * m * atDeltas), atDeltas)
    return (qn0, qn1)

def BetaMN(qNME, qNMI, n, m):
    # qMNE and qMNI are of size 2-by-M
    qnm = np.array([qNME[n, m], qNMI[n, m]])
    return np.dot(Jab **2, qnm)

def GenQuenchedNoiseWithCorr(n):
    x = np.random.randn(n)
    y = np.random.randn(n)
    # print choleskyMatE
    v = np.dot(choleskyMatE[:n, :n], x)
    w = np.dot(choleskyMatE[:n, :n], y)
    out = lambda atPhi: np.dot(v, np.cos(2 * np.arange(n) * atPhi)) + np.dot(w, np.sin(2 * np.arange(n) * atPhi))
    return out

def OSI(firingRate, atTheta):
    # theta in radians
    out = np.nan
    zk = np.dot(firingRate, np.exp(2j * atTheta))
    if(firingRate.mean() > 0.0):
        out = np.absolute(zk) / np.sum(firingRate)
    return out

def OSIOfPop(firingRates, atThetas):
    # thetas in radians
    nNeurons, nThetas = firingRates.shape
    out = np.zeros((nNeurons, ))
    for i in range(nNeurons):
        out[i] = OSI(firingRates[i , :], atThetas)
    return out

if __name__ == "__main__":
    #[jee, Jei
    # jie, jii]
    Jab = np.array([[1.0, -1.5],
                    [1.0, -1.00]])
    m0 = 0.075
    mExtZero = m0
    mExtOne = 0.075
    cFF = 0.20
    Ea = np.array([2.0, 1.0]) #/ cFF
    p = np.array([0, 0])
    gamma = np.array([0, 0]) #  gammaE and gammaI 
    mFourier0 = -1.0 * np.dot(np.linalg.inv(Jab), Ea) * m0
    mFourier1 = np.array([0.05, 0.01])

    uInitialGuess = [0.1, 0.1] 
    u0Solution = fsolve(meanInput, uInitialGuess, args = (mFourier0[0], mFourier1[0]), full_output = 1)
    PrintError(u0Solution, 'uE')
    u1Solution = fsolve(meanInput, uInitialGuess, args = (mFourier0[1], mFourier1[1]), full_output = 1)
    PrintError(u1Solution, 'uI')
    u0E, u1E = u0Solution[0]
    u0I, u1I = u1Solution[0]
    # outQ = fsolve(QuenchedDisorder, (0.01, 0.03478, 0.00, 0.00), args = (mFourier0, mFourier1, u0E, u0I, u1E, u1I), xtol = 1e-10, full_output = 1)
    
    # PrintError(outQ, 'q0 q1')
    
    # q0E, q0I, q1E, q1I =  fsolve(QuenchedDisorder, (0.01, 0.02, 0.00, 0.00), args = (mFourier1, mFourier1, u0E, u0I, u1E, u1I, 0))


    q0E, q0I, q1E, q1I =  fsolve(QuenchedDisorder, (0.02, 0.05, 0.001, 0.001), args = (mFourier1, mFourier1, u0E, u0I, u1E, u1I, 0))
    
    # q0E, q0I, q1E, q1I = outQ[0]

    betaAtPhi = lambda atPhi: np.dot(Jab **2, np.array([q0E, q0I]) + p[0] * np.array([q1E, q1I]) * np.cos(2 * atPhi))
    alphaA = np.dot(Jab**2, mFourier0 + p * mFourier1 * np.cos(2 * 0))

    print "\n\nSummary:"   
    print "-------------------------------"
    print "           E         |    I"
    print "-------------------------------"
    print "m0    : %.2e       %.2e"%(mFourier0[0], mFourier0[1])
    print "m1    : %.2e       %.2e"%(mFourier1[0], mFourier1[1])    
    print "u0    : %.2e       %.2e"%(u0E, u0I)
    print "u1    : %.2e       %.2e"%(u1E, u1I)
    print "q0    : %.2e       %.2e"%(q0E, q0I)
    print "q1    : %.2e       %.2e"%(q1E, q1I)
    print "alpha : %.2e       %.2e"%(alphaA[0], alphaA[1])
    print "beta  : %.2e       %.2e"%(betaAtPhi(0)[0], betaAtPhi(0)[1])    

    if betaAtPhi(0)[0] > alphaA[0]:
	print 'beta > alpha!!!!!'
	raise SystemExit

    nNeurons = 100
    atDeltas = np.linspace(0, np.pi, 10)
    atPhis = np.linspace(-np.pi/2, np.pi/2, 17)
    IF_COMPUTE = 1
    maxCorrM = 5
    if IF_COMPUTE:
        qnOfDeltaE = np.zeros((2, atDeltas.size))
        qnOfDeltaI = np.zeros((2, atDeltas.size))
        # bbplt = np.zeros((atDeltas.size, ))
        for k, kDelta in enumerate(atDeltas):
            tmp = FSolveAtDelta2(kDelta)
            qq0 = tmp[0]
            qq1 = tmp[1]
            qnOfDeltaE[0, k] = qq0[0]
            qnOfDeltaE[1, k] = qq1[0]
            qnOfDeltaI[0, k] = qq0[1]
            qnOfDeltaI[1, k] = qq1[1]
            # bbplt[k] = np.dot(Jab[0, :] **2, np.array(qq0) + p * np.array(qq1))
            maxM = 25
            qnmE = np.zeros((2, maxM))
            qnmI = np.zeros((2, maxM))
	    # ipdb.set_trace()
            for mm in range(maxM):
                qnmE[:, mm] = QMN(mm, qnOfDeltaE, 0)
                qnmI[:, mm] = QMN(mm, qnOfDeltaI, 1)
            betaNME = np.zeros((2, maxM))
            betaNMI = np.zeros((2, maxM))
            for nn in range(2):
                for mm in range(maxM):
                    betaNME[nn, mm]  = BetaMN(qnmE, qnmI, nn, mm)[0]
                    betaNMI[nn, mm]  = BetaMN(qnmE, qnmI, nn, mm)[1]
        tmpQQQxxx0 = np.array([qnOfDeltaE[0, 0], qnOfDeltaI[0, 0]])
        tmpQQQxxx1 = np.array([qnOfDeltaE[1, 0], qnOfDeltaI[1, 0]])
        print "xxxxxxxxxxxxxxxxxxxxxxx"
        print tmpQQQxxx0
        print tmpQQQxxx1
        print np.sum(np.sum(betaNME)), np.dot(np.dot(np.array([1., 1.]), betaNME), np.cos(2 * 0 * np.arange(maxM))), np.dot(Jab[0, :]**2, (tmpQQQxxx0 + p * tmpQQQxxx1))
        print "--------------------------"
        qCorrMatE = np.zeros((maxCorrM, maxCorrM))
        qCorrMatI = np.zeros((maxCorrM, maxCorrM))
        for n in range(maxCorrM):
            for m in range(maxCorrM):
                if(np.abs(n - m) < 2 and (n + m) < maxM):
                    qCorrMatE[n, m] = betaNME[np.abs(n - m), n + m]
        np.save('corrmat', qCorrMatE)
    qCorrMatE = np.load('corrmat.npy')
    print "EIGS--->"
    W, V = np.linalg.eig(qCorrMatE)
    print W
    choleskyMatE =  np.linalg.cholesky(qCorrMatE)

    # _, osi = TuningCurve(nNeurons = 500000)
    # print "osi analytic: ", np.nanmean(osi)
    # plt.hist(osi, 100, normed = 1, histtype = 'step', color = 'c')

    
    # sc.PltOSIHist(p[0], gamma[0], 8, mExtZero, mExtOne, IF_NEW_FIG = False, K = 500, T=2000, color = 'r')
    sc.PltOSIHist(p[0], gamma[0], 8, mExtZero, mExtOne, IF_NEW_FIG = False, color = 'g')
    # sc.PltOSIHist(p[0], gamma[0], 8, mExtZero, mExtOne, IF_NEW_FIG = False, K = 1500, T=2000, color = 'b')
    sc.PltOSIHist(0, 0, 36, mExtZero, mExtOne, K=2000, T=10000, trNo=11, IF_NEW_FIG = False, color = 'k')
    # plt.legend(['analytic', 'K=500', 'K=1000', 'K=1500', 'K=2000'], loc = 0, frameon = False, numpoints = 1, prop = {'size': 10})
    plt.legend(['analytic', 'K=1000', 'K=2000'], loc = 0, frameon = False, numpoints = 1, prop = {'size': 10})    
    plt.savefig('./figs/twopop/osi_distr_theory_sim.png')
    # plt.show()

   #####
    kOSI = np.zeros((10, ))
    pool = Pool(processes=12)
    for k in range(1000):
    	print k, ' ',
    	sys.stdout.flush()
    	kOSI[k] = np.nanmean(TuningCurve(nNeurons = 10000))

    # result = pool.apply_async(TuningCurve, [100])
    # ipdb.set_trace()
    # kOSI = np.array(result)
    # kOSI = np.array(kOSI)
    kSamples = np.sum(~np.isnan(kOSI))
    print 'osi analytic = ', np.nanmean(kOSI), ', SEM = ', np.nanstd(kOSI) / np.sqrt(kSamples)
    plt.figure()
    np.save('mean_osi_distr_samples_new', kOSI)
    plt.hist(kOSI, 2, normed = 1, histtype = 'step', color = 'k')
    plt.savefig('./figs/twopop/osi_mean_distr_n1000.png')
    plt.show()

