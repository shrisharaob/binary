from __future__ import division
import numpy as np
# import mpmath as mp
import pylab as plt
from scipy.optimize import fsolve, brenth, root
from scipy.integrate import quad, simps
from scipy.special import erfc 
from multiprocessing import Pool
from functools import partial 
import sys

def MyErfc(z):
    # normal CDF
    return 0.5 * erfc(z / np.sqrt(2.0))

def MyErfc3(z):
    # normal CDF sqrd
    return (0.5 * erfc(z / np.sqrt(2.0))) **2


def MyErfc2(z):
    # normal CDF sqrd
    return (MyErfc(z)) **2        

def alphaF(phi, EorI):
    out = np.dot(Jab**2, mFourier0 + p * mFourier1 * np.cos(2 * phi))    
    return np.abs(out[EorI])

def betaOfPhi(phi, delta, q0, q1, EorI):
    out = np.dot(Jab **2, q0 + p * q1 * np.cos(2 * (phi + delta)))    
    return np.abs(out[EorI])

def meanInput(uintial, *args):
    u0, u1 = uintial
    mFourier0, mFourier1 = args
    func0 = lambda phi: ((1.0 / (np.pi)) * MyErfc((1.0 - u0 - u1 * np.cos(2 * phi)) / np.sqrt(alphaF(phi, EorI = 0))))
    func1 = lambda phi: ((2.0 / (np.pi)) * MyErfc((1.0 - u0 - u1 * np.cos(2 * phi)) / np.sqrt(alphaF(phi, EorI = 1))) * np.cos(2 * phi))
    a = mFourier0 - quad(func0, 0, np.pi)[0]
    b = mFourier1 - quad(func1, 0, np.pi)[0]
    return (a, b)

def CheckQ(q0):
    if q0[0] > mFourier0[0]:
        q0[0] = mFourier0[0] - 1e-2
    if q0[1] > mFourier0[1]:
        q0[1] = mFourier0[1] - 1e-2
    # print q0        
    return q0
        
    
def quenchedDisorder(qinitial, *args):
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
                                # func1E = lambda phi: (2.0 / np.pi) * np.dot(MyErfc((1 - uFuncE(phi, 0) - np.sqrt(betaOfPhi(phi, 0, q0, q1, EorI = 0)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, 0, q0, q1, EorI = 0))) * MyErfc((1 - uFuncE(phi, 0) - np.sqrt(betaOfPhi(phi, 0, q0, q1, EorI = 0)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, 0, q0, q1, EorI = 0))) * np.cos(2 * phi), hermiteWeights)     
    
    func1I = lambda phi: (2.0 / np.pi) * np.dot(MyErfc((1 - uFuncI(phi, delta) - np.sqrt(betaOfPhi(phi, delta, q0, q1, EorI = 1)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, delta, q0, q1, EorI = 1))) * MyErfc((1 - uFuncI(phi, -delta) - np.sqrt(betaOfPhi(phi, -delta, q0, q1, EorI = 1)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, -delta, q0, q1, EorI = 1)))  * np.cos(2 * phi), hermiteWeights) #* np.cos(2 * phi)
    
    
    # func0E = lambda phi: (1.0 / (np.pi * np.sqrt(np.pi))) * np.dot(erfc((1 - uFuncE(phi) - np.sqrt(betaOfPhi(phi, q0, q1, EorI = 0)) * hermiteX * np.sqrt(2.0)) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, q0, q1, EorI = 0))) **2, hermiteWeights)
    # func0I = lambda phi: (1.0 / (np.pi * np.sqrt(np.pi))) * np.dot(erfc((1 - uFuncI(phi) - np.sqrt(betaOfPhi(phi, q0, q1, EorI = 1)) * hermiteX * np.sqrt(2.0)) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, q0, q1, EorI = 1))) **2, hermiteWeights)

    # func1E = lambda phi: (2.0 / (np.pi * np.sqrt(np.pi))) * np.dot(erfc((1 - uFuncE(phi) - np.sqrt(betaOfPhi(phi, q0, q1, EorI = 0)) * hermiteX * np.sqrt(2.0)) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, q0, q1, EorI = 0))) **2, hermiteWeights) * np.cos(2 * phi)
    # func1I = lambda phi: (2.0 / (np.pi * np.sqrt(np.pi))) * np.dot(erfc((1 - uFuncI(phi) - np.sqrt(betaOfPhi(phi, q0, q1, EorI = 1)) * hermiteX * np.sqrt(2.0)) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, q0, q1, EorI = 1))) **2, hermiteWeights) * np.cos(2 * phi)
    
    # func1E = lambda phi: (2.0 / (np.pi * np.sqrt(np.pi))) * np.dot((erfc((1 - uFuncE(phi) - np.sqrt(betaOfPhi(phi, q0, q1, EorI = 0)) * hermiteX * np.sqrt(2.0)) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, q0, q1, EorI = 0))) **2) * np.cos(2 * phi), hermiteWeights)
    # func1I = lambda phi: (2.0 / (np.pi * np.sqrt(np.pi))) * np.dot((erfc((1 - uFuncI(phi) - np.sqrt(betaOfPhi(phi, q0, q1, EorI = 1)) * hermiteX * np.sqrt(2.0)) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, q0, q1, EorI = 1))) **2)  * np.cos(2 * phi) , hermiteWeights) 

    # aE = quad(func0E, 0, np.pi)[0] - np.abs(q0E)
    # aI = quad(func0I, 0, np.pi)[0] - np.abs(q0I)
    # bE = quad(func1E, 0, np.pi)[0] - np.abs(q1E)
    # bI = quad(func1I, 0, np.pi)[0] - np.abs(q1I)

    # print "::", alphaF(0, EorI = 0), betaOfPhi(0, 0, q0, q1, EorI = 1), quad(func1E, 0, np.pi)[0],  q1E
    
    aE = quad(func0E, 0, np.pi)[0] - q0E
    aI = quad(func0I, 0, np.pi)[0] - q0I
    bE = quad(func1E, 0, np.pi)[0] - q1E
    bI = quad(func1I, 0, np.pi)[0] - q1I
    return (aE, aI, bE, bI)


def FSolveAvgAtPhi(atPhi):
    print "phi = ", atPhi * 180. / np.pi
    alphaAtPhi = np.dot(Jab**2, mFourier0 + p * mFourier1 * np.cos(2 * atPhi))    
    u0E, u1E =  fsolve(meanInput, uInitialGuess, args = (mFourier0[0], mFourier1[0]))
    u0I, u1I =  fsolve(meanInput, uInitialGuess, args = (mFourier0[1], mFourier1[1]))
    return erfc((1.0 - u0E - u1E * np.cos(2 * atPhi)) / np.sqrt(alphaAtPhi[0]))

def FSolveAtPhi(atPhi):
    print "phi = ", atPhi * 180. / np.pi
    alphaAtPhi = np.dot(Jab**2, mFourier0 + p * mFourier1 * np.cos(2 * atPhi))
    u0E, u1E =  fsolve(meanInput, uInitialGuess, args = (mFourier0[0], mFourier1[0]))
    u0I, u1I =  fsolve(meanInput, uInitialGuess, args = (mFourier0[1], mFourier1[1]))    
    q0E, q0I, q1E, q1I =  fsolve(quenchedDisorder, (0.010974952807053, 0.03477881883636, 0.00, 0.00), args = (mFourier0, mFourier1, u0E, u0I, u1E, u1I), xtol = 1e-10)
    q0 = np.array([q0E, q0I])
    q1 = np.array([q1E, q1I])
    betaAtPhi = np.dot(Jab **2, q0 + p * q1 * np.cos(2 * atPhi))
    nNeurons = 1000
    mi = np.zeros((nNeurons, ))
    nIterations = 1
#    for iii in range(nIterations):
    xi = np.sort(np.random.randn(nNeurons))
    for i in np.arange(nNeurons):
        mi[i] = erfc((1.0 - u0E - u1E * np.cos(2 * atPhi) - np.sqrt(betaAtPhi[0]) * xi[i]) / np.sqrt(alphaAtPhi[0] - betaAtPhi[0]))
    mOfPhi = np.mean(mi)
    return mOfPhi #(u0E, u0I, u1E, u1I, q0E, q0I, q1E, q1I, alphaAtPhi, betaAtPhi)

def FSolveAtDelta(atPhi, atDelta):
    # print "delta = ", atDelta * 180. / np.pi
   # alphaAtPhi = np.dot(Jab**2, mFourier0 + p * mFourier1 * np.cos(2 * atPhi))
    u0E, u1E =  fsolve(meanInput, uInitialGuess, args = (mFourier0[0], mFourier1[0]))
    u0I, u1I =  fsolve(meanInput, uInitialGuess, args = (mFourier0[1], mFourier1[1]))        
  #  atDeltas = np.linspace(0, np.pi/2, 4)
    q0E, q0I, q1E, q1I =  fsolve(quenchedDisorder, (0.01, 0.02, 0.00, 0.00), args = (mFourier1, mFourier1, u0E, u0I, u1E, u1I, atDelta))
    q0 = np.array([q0E, q0I])
    q1 = np.array([q1E, q1I])    
    betaAtDelta = np.dot(Jab **2, q0 + p * q1 * np.cos(2 * atPhi))
    return  betaAtDelta #(u0E, u0I, u1E, u1I, q0E, q0I, q1E, q1I, alphaAtPhi, betaAtPhi)

def FSolveAtDelta2(atDelta):
    # Returns q(0)OfDelta and q(1)OfDelta
    print atDelta * 180.000 / np.pi
    u0E, u1E =  fsolve(meanInput, uInitialGuess, args = (mFourier0[0], mFourier1[0]))
    u0I, u1I =  fsolve(meanInput, uInitialGuess, args = (mFourier0[1], mFourier1[1]))

    # outQ = fsolve(quenchedDisorder, (0., 0., 0.00, 0.00), args = (mFourier0, mFourier1, u0E, u0I, u1E, u1I, 0), xtol = 1e-10, full_output = 1)
    # PrintError(outQ, 'q0 q1')
    
    q0E, q0I, q1E, q1I =  fsolve(quenchedDisorder, (0.01, 0.02, 0.00, 0.00), args = (mFourier1, mFourier1, u0E, u0I, u1E, u1I, atDelta))
    q0 = np.array([q0E, q0I])
    q1 = np.array([q1E, q1I])
    # print q0, q1
    return  (q0, q1)

def GenQuenchedNoiseWithCorr(n):
    x = np.random.randn(n)
    y = np.random.randn(n)
    # print choleskyMatE
    v = np.dot(choleskyMatE[:n, :n], x)
    w = np.dot(choleskyMatE[:n, :n], y)
    out = lambda atPhi: np.dot(v, np.cos(2 * np.arange(n) * atPhi)) + np.dot(w, np.sin(2 * np.arange(n) * atPhi))
    return out
    
    
def TuningCurve():
    u0E, u1E =  fsolve(meanInput, uInitialGuess, args = (mFourier0[0], mFourier1[0]))
    u0I, u1I =  fsolve(meanInput, uInitialGuess, args = (mFourier0[1], mFourier1[1]))
    q0E, q0I, q1E, q1I =  fsolve(quenchedDisorder, (0.01, 0.02, 0.00, 0.00), args = (mFourier1, mFourier1, u0E, u0I, u1E, u1I, 0))
    q0 = np.array([q0E, q0I])
    q1 = np.array([q1E, q1I])    
    betaAtDeltaZero = lambda atPhi: betaOfPhi(atPhi, 0, q0, q1, EorI = 0)
    # atPhis = np.linspace(-np.pi/2, np.pi/2, 9)
    miOfPhi = np.zeros((atPhis.size, ))
    # pool = Pool(atPhis.size)
    nNeurons = 100
    plt.ion()

    for l in range(nNeurons):
        qNoise = GenQuenchedNoiseWithCorr(maxCorrM)        
        for i, iPhi in enumerate(atPhis):
            miOfPhi[i] = MyErfc((1.0 - u0E - u1E * np.cos(2 * iPhi) - qNoise(iPhi)) / np.sqrt(alphaF(iPhi, 0) - betaAtDeltaZero(iPhi)))
            # miOfPhi[i] = MyErfc((1.0 - u0E - u1E * np.cos(2 * iPhi)) / np.sqrt(alphaF(iPhi, 0) - betaAtDeltaZero(iPhi)))
        # print miOfPhi
        plt.plot(atPhis * 180 / np.pi, miOfPhi, 'ko-')
        plt.grid()
        plt.xlabel(r'$\phi$(deg)', fontsize = 20)
        plt.ylabel(r'$m^i_E(\phi)$', fontsize = 20)
        plt.savefig('fig%s'%(l))
        # plt.show()
        plt.waitforbuttonpress()
        plt.clf()
    return miOfPhi
    
def FourierOfBeta(beta, phi, delta, n, m):
    #beta has to be function handle
    lenPhiSummation = phi.size
    lenDeltaSummation = delta.size
    dPhi = phi[1] - phi[0]
    dDelta = delta[1] - delta[0]
    out = 0;
    for i in range(lenPhiSummation):
        for j in range(lenDeltaSummation):
            out += dPhi * dDelta * beta[i, j] * np.cos(2 * n * phi[i]) * np.cos(2 * m * delta[j])
    return out

def FourierOfBeta2(beta, phi, delta, n, m):
    #beta has to be function handle
    lenPhiSummation = phi.size
    lenDeltaSummation = delta.size
    dPhi = phi[1] - phi[0]
    dDelta = delta[1] - delta[0]
    out = 0;
    for i in range(lenPhiSummation):
        for j in range(lenDeltaSummation):
            out += dPhi * dDelta * beta[i, j] * np.cos(2 * n * phi[i]) * np.cos(2 * m * delta[j])
    return out

def QMN(m, qnOfDelta):
    # qnOfDelta must be 2-by-nDelta matrix
    if m == 0:
        qn0 = (1. / np.pi) * simps(qnOfDelta[0, :], atDeltas)
        qn1 = (p / np.pi) * simps(qnOfDelta[1, :], atDeltas)        
    else:
        qn0 = (2. / np.pi) * simps(qnOfDelta[0, :] * np.cos(2 * m * atDeltas), atDeltas)
        qn1 = (p * 2.0 / np.pi) * simps(qnOfDelta[1, :] * np.cos(2 * m * atDeltas), atDeltas)
    # print qm0, qm1
    return (qn0, qn1)

def BetaMN(qNME, qNMI, n, m):
    # qMNE and qMNI are of size 2-by-M
    qnm = np.array([qNME[n, m], qNMI[n, m]])
    return np.dot(Jab **2, qnm)

def PrintError(outStruct, varName):
    print "-----------------------------------"
    print "Variable        :",  varName
    print "Solution        :", outStruct[0]
    print "error           :", outStruct[1]['fvec']
    print "Solution Status :", outStruct[-1]
    print "-----------------------------------"    


if __name__ == "__main__":
    Jab = np.array([[1.0, -1.5],
                    [1.0, -1.00]])
    m0 = 0.075
    Ea = np.array([2.0, 1.0])
    p = 0.25 # MAX is 0.5 
    mu = 0.25
    mFourier0 = -1.0 * np.dot(np.linalg.inv(Jab), Ea) * m0
    if p == 0:
        mFourier1 = np.array([0, 0])
    else:
        mFourier1 = -1.0 *  (mu / p) * np.dot(np.linalg.inv(Jab), Ea) * m0

    uInitialGuess = [0.1, 0.1] 
    # u0Solution = fsolve(meanInput, uInitialGuess, args = (mFourier0[0], mFourier1[0]), full_output = 1)
    # PrintError(u0Solution, 'uE')
    # u1Solution = fsolve(meanInput, uInitialGuess, args = (mFourier0[1], mFourier1[1]), full_output = 1)
    # PrintError(u1Solution, 'uI')
    # u0E, u1E = u0Solution[0]
    # u0I, u1I = u1Solution[0]

    # outQ = fsolve(quenchedDisorder, (0.010974952807053, 0.03477881883636, 0.00, 0.00), args = (mFourier0, mFourier1, u0E, u0I, u1E, u1I), xtol = 1e-10, full_output = 1)
    # PrintError(outQ, 'q0 q1')
    # q0E, q0I, q1E, q1I = outQ[0]

    # betaAtPhi = lambda atPhi: np.dot(Jab **2, np.array([q0E, q0I]) + p * np.array([q1E, q1I]) * np.cos(2 * atPhi))
    alphaA = np.dot(Jab**2, mFourier0 + p * mFourier1 * np.cos(2 * 0))
    
    # print "\n\nSummary:"   
    # print "-------------------------------"
    # print "           E         |    I"
    # print "-------------------------------"
    # print "m0    : %.2e       %.2e"%(mFourier0[0], mFourier0[1])
    # print "m1    : %.2e       %.2e"%(mFourier1[0], mFourier1[1])    
    # print "u0    : %.2e       %.2e"%(u0E, u0I)
    # print "u1    : %.2e       %.2e"%(u1E, u1I)
    # print "q0    : %.2e       %.2e"%(q0E, q0I)
    # print "q1    : %.2e       %.2e"%(q1E, q1I)
    # print "alpha : %.2e       %.2e"%(alphaA[0], alphaA[1])
    # print "beta  : %.2e       %.2e"%(betaAtPhi(0)[0], betaAtPhi(0)[1])    


    
    # print "\n\nalpha > beta:", np.all(alphaA > betaAtPhi(0))


    # PLOT
    IF_PLOT = 0
    if len(sys.argv) > 1:
        IF_PLOT = int(sys.argv[1])
    if IF_PLOT:
        atPhis = np.linspace(-np.pi/2, np.pi/2, 17)
        mOfPhi = np.zeros((atPhis.size, ))
        avgmOfPhi = np.zeros((atPhis.size, ))
        pool = Pool(atPhis.size)
        avgmOfPhi = pool.map(FSolveAtPhi, atPhis)
        mOfPhi0 = np.zeros((atPhis.size, ))
        mOfPhi0 = np.array(pool.map(FSolveAvgAtPhi, atPhis))
        pool.close()

        plt.plot(atPhis * 180.0 / np.pi, avgmOfPhi, 'ko-')
        plt.plot(atPhis * 180.0 / np.pi, mOfPhi0, 'ro-')
        plt.xlabel(r'$\phi$(deg)', fontsize = 20)
        plt.ylabel(r'$m_E(\phi)$', fontsize = 20)
        plt.show()

    nNeurons = 100
    atDeltas = np.linspace(0, np.pi, 129)
    # atDeltas = np.array([0, np.pi])
    atPhis = np.linspace(-np.pi/2, np.pi/2, 17)
    # atPhis = np.array([0, np.pi / 2, np.pi])

    # print atDeltas.size
    IF_COMPUTE = 1
    maxCorrM = 5
    if IF_COMPUTE:
        # pool = Pool(int(atDeltas.size))
        # betaMatE = np.zeros((atPhis.size, atDeltas.size))
        # betaMatI = np.zeros((atPhis.size, atDeltas.size))
        # betaMNE = np.zeros((maxN, maxM))
        # betaMNI = np.zeros((maxN, maxM))

        qnOfDeltaE = np.zeros((2, atDeltas.size))
        qnOfDeltaI = np.zeros((2, atDeltas.size))


        bbplt = np.zeros((atDeltas.size, ))
        
        for k, kDelta in enumerate(atDeltas):
            tmp = FSolveAtDelta2(kDelta)
            qq0 = tmp[0]
            qq1 = tmp[1]
            qnOfDeltaE[0, k] = qq0[0]
            qnOfDeltaE[1, k] = qq1[0]
            qnOfDeltaI[0, k] = qq0[1]
            qnOfDeltaI[1, k] = qq1[1]
            bbplt[k] = np.dot(Jab[0, :] **2, np.array(qq0) + p * np.array(qq1))


            # plt.plot(atDeltas, bbplt)


            # print "q(0)", qnOfDeltaE[:, 0]
            # print qnOfDeltaI[:, 0]            

            maxM = 25
            qnmE = np.zeros((2, maxM))
            qnmI = np.zeros((2, maxM))
            for mm in range(maxM):
                qnmE[:, mm] = QMN(mm, qnOfDeltaE)
                qnmI[:, mm] = QMN(mm, qnOfDeltaI)

            betaNME = np.zeros((2, maxM))
            betaNMI = np.zeros((2, maxM))
            for nn in range(2):
                for mm in range(maxM):
                    betaNME[nn, mm]  = BetaMN(qnmE, qnmI, nn, mm)[0]
                    betaNMI[nn, mm]  = BetaMN(qnmE, qnmI, nn, mm)[1]

        # zzz = np.sum(betaNME, 0)
        # yyy = lambda dlta: np.sum(zzz * np.cos(2 * range(10) * dlta)
        tmpQQQxxx0 = np.array([qnOfDeltaE[0, 0], qnOfDeltaI[0, 0]])
        tmpQQQxxx1 = np.array([qnOfDeltaE[1, 0], qnOfDeltaI[1, 0]])

        print "xxxxxxxxxxxxxxxxxxxxxxx"
        print tmpQQQxxx0
        print tmpQQQxxx1
        print " YOYOYO"
        print np.sum(np.sum(betaNME)), np.dot(np.dot(np.array([1., 1.]), betaNME), np.cos(2 * 0 * np.arange(maxM))), np.dot(Jab[0, :]**2, (tmpQQQxxx0 + p * tmpQQQxxx1))



        print "--------------------------"

        qCorrMatE = np.zeros((maxCorrM, maxCorrM))
        qCorrMatI = np.zeros((maxCorrM, maxCorrM))
        for n in range(maxCorrM):
            for m in range(maxCorrM):
                if(np.abs(n - m) < 2 and (n + m) < maxM):
                    qCorrMatE[n, m] = betaNME[np.abs(n - m), n + m]
                    # qCorrMatE[n, m] = betaNME[np.abs(n - m), n + m]                    

        np.save('corrmat', qCorrMatE)

    qCorrMatE = np.load('corrmat.npy')

    print "EIGS--->"
    W, V = np.linalg.eig(qCorrMatE)
    print W

    choleskyMatE =  np.linalg.cholesky(qCorrMatE)
    # plt.matshow(choleskyMatE[:4, :4])
    # print choleskyMatE[:4, :4]
    # plt.show()

    TuningCurve()
                
            
            
        
        



            

        
        # # atPhis = np.array([0])
        # for ii, iPhi in enumerate(atPhis):
        #     print ii
        #     funcFun = partial(FSolveAtDelta, iPhi)
        #     tmp = np.asarray(pool.map(funcFun, atDeltas))
        #     qq0 = tmp[0]
        #     qq1 = tmp[1]
        #     qnOfDeltaE[:, 
        #     # print tmp
        #     betaMatE[ii, :] = tmp[:, 0]
        #     betaMatI[ii, :] = tmp[:, 1]
        # print betaMatE
        # pool.close()            
        # for nn in range(maxN):
        #     for mm in range(maxM):
        #         betaMNE[nn, mm] = FourierOfBeta(betaMatE, atPhis, atDeltas, nn, mm)
        #         betaMNI[nn, mm] = FourierOfBeta(betaMatI, atPhis, atDeltas, nn, mm)            
        # np.save('betas', [betaMatE, betaMatI, betaMNE, betaMNI])
        # print betaMNE
    # betaMatE, betaMatI, betaMNE, betaMNI = np.load('betas.npy')
    # print betaMNE.shape
    
    # qMaxN = 2
    # qCorrMatE = np.zeros((qMaxN, qMaxN))
    # qCorrMatI = np.zeros((qMaxN, qMaxN))    
    # for n in range(maxN):
    #     for m in range(qMaxN):
    #         qCorrMatE[n, m] = betaMNE[np.abs(n - m), n + m]

    # choleskyMatE =  np.linalg.cholesky(qCorrMatE)



    # TuningCurve()
            
    
    

      
        

