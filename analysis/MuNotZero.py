from __future__ import division
import numpy as np
import mpmath as mp
import pylab as plt
from scipy.optimize import fsolve, brenth, root
from scipy.integrate import quad
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

def betaOfPhi2(phi, q0, q1):
    out = np.dot(Jab **2, q0 + p * q1 * np.cos(2 * phi))
    print out
    return np.abs(out)
    
def betaOfPhi(phi, q0, q1, EorI):
    out = np.dot(Jab **2, q0 + p * q1 * np.cos(2 * phi))    
    # localAlpha = np.dot(Jab**2, mFourier0 + p * mFourier1 * np.cos(2 * phi))[EorI]
    # if localAlpha <= out[EorI]:
    #     out[EorI] = localAlpha - 1e-6
    #print out
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
    mFourier0, mFourier1, u0E, u0I, u1E, u1I = args
    q0 = np.array([q0E, q0I])
    q1 = np.array([q1E, q1I])
    hermiteDeg = 30
    hermiteX, hermiteWeights = np.polynomial.hermite.hermgauss(hermiteDeg)
    hermiteX = hermiteX * np.sqrt(2.0)
    hermiteWeights = hermiteWeights / np.sqrt(np.pi)
    uFuncE = lambda phi: (u0E + u1E * np.cos(2 * phi))
    uFuncI = lambda phi: (u0I + u1I * np.cos(2 * phi))    
    func0E = lambda phi: (1.0 / np.pi) * np.dot(MyErfc2((1 - uFuncE(phi) - np.sqrt(betaOfPhi(phi, q0, q1, EorI = 0)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, q0, q1, EorI = 0))), hermiteWeights)
    func0I = lambda phi: (1.0 / np.pi) * np.dot(MyErfc2((1 - uFuncI(phi) - np.sqrt(betaOfPhi(phi, q0, q1, EorI = 1)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, q0, q1, EorI = 1))), hermiteWeights)
    func1E = lambda phi: (2.0 / np.pi) * np.dot(MyErfc2((1 - uFuncE(phi) - np.sqrt(betaOfPhi(phi, q0, q1, EorI = 0)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 0) - betaOfPhi(phi, q0, q1, EorI = 0))) * np.cos(2 * phi), hermiteWeights) 
    func1I = lambda phi: (2.0 / np.pi) * np.dot(MyErfc2((1 - uFuncI(phi) - np.sqrt(betaOfPhi(phi, q0, q1, EorI = 1)) * hermiteX) / np.sqrt(alphaF(phi, EorI = 1) - betaOfPhi(phi, q0, q1, EorI = 1))) * np.cos(2 * phi), hermiteWeights) 
    aE = quad(func0E, 0, np.pi)[0] - q0E
    aI = quad(func0I, 0, np.pi)[0] - q0I
    bE = quad(func1E, 0, np.pi)[0] - q1E
    bI = quad(func1I, 0, np.pi)[0] - q1I
    return (aE, aI, bE, bI)

def FSolveAvgAtPhi(atPhi):
 #   print "phi = ", atPhi * 180. / np.pi
    alphaAtPhi = np.dot(Jab**2, mFourier0 + p * mFourier1 * np.cos(2 * atPhi))    
    u0E, u1E =  fsolve(meanInput, uInitialGuess, args = (mFourier0[0], mFourier1[0]))
    u0I, u1I =  fsolve(meanInput, uInitialGuess, args = (mFourier0[1], mFourier1[1]))
    return MyErfc((1.0 - u0E - u1E * np.cos(2 * atPhi)) / np.sqrt(alphaAtPhi[0]))

def FSolveAtPhi(atPhi):
#    print "phi = ", atPhi * 180. / np.pi
    alphaAtPhi = np.dot(Jab**2, mFourier0 + p * mFourier1 * np.cos(2 * atPhi))
    u0E, u1E =  fsolve(meanInput, uInitialGuess, args = (mFourier0[0], mFourier1[0]))
    u0I, u1I =  fsolve(meanInput, uInitialGuess, args = (mFourier0[1], mFourier1[1]))    
    q0E, q0I, q1E, q1I =  fsolve(quenchedDisorder, (0.010974952807053, 0.03477881883636, 0.00, 0.00), args = (mFourier0, mFourier1, u0E, u0I, u1E, u1I), xtol = 1e-10)
    q0 = np.array([q0E, q0I])
    q1 = np.array([q1E, q1I])
    betaAtPhi = np.dot(Jab **2, q0 + p * q1 * np.cos(2 * atPhi))
    nNeurons = 50000
    mi = np.zeros((nNeurons, ))
    xi = np.sort(np.random.randn(nNeurons))
    for i in np.arange(nNeurons):
        mi[i] = MyErfc((1.0 - u0E - u1E * np.cos(2 * atPhi) - np.sqrt(betaAtPhi[0]) * xi[i]) / np.sqrt(alphaAtPhi[0] - betaAtPhi[0]))
    mOfPhi = np.mean(mi)
    return mOfPhi #(u0E, u0I, u1E, u1I, q0E, q0I, q1E, q1I, alphaAtPhi, betaAtPhi)


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
    p = 0.25
    mu = 0.05
    mFourier0 = -1.0 * np.dot(np.linalg.inv(Jab), Ea) * m0
    if p == 0:
        mFourier1 = np.array([0, 0])
    else:
        mFourier1 = -1.0 *  (mu / p) * np.dot(np.linalg.inv(Jab), Ea) * m0

    uInitialGuess = [0.1, 0.1] 
    u0Solution = fsolve(meanInput, uInitialGuess, args = (mFourier0[0], mFourier1[0]), full_output = 1)
    PrintError(u0Solution, 'uE')
    u1Solution = fsolve(meanInput, uInitialGuess, args = (mFourier0[1], mFourier1[1]), full_output = 1)
    PrintError(u1Solution, 'uI')
    u0E, u1E = u0Solution[0]
    u0I, u1I = u1Solution[0]

#    outQ = fsolve(quenchedDisorder, (0.010974952807053, 0.03477881883636, 0.00, 0.00), args = (mFourier0, mFourier1, u0E, u0I, u1E, u1I), xtol = 1e-10, full_output = 1)
    outQ = fsolve(quenchedDisorder, (0, 0, 0.00, 0.00), args = (mFourier0, mFourier1, u0E, u0I, u1E, u1I), xtol = 1e-10, full_output = 1)    
    PrintError(outQ, 'q0 q1')
    q0E, q0I, q1E, q1I = outQ[0]

    betaAtPhi = lambda atPhi: np.dot(Jab **2, np.array([q0E, q0I]) + p * np.array([q1E, q1I]) * np.cos(2 * atPhi))
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


    
    print "\n\nalpha > beta:", np.all(alphaA > betaAtPhi(0))


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

    plt.plot(atPhis * 180.0 / np.pi, avgmOfPhi, 'ks-', markersize = 6, label = 'mean')
    plt.plot(atPhis * 180.0 / np.pi, mOfPhi0, 'ro-', label = 'avgOverRealizations')
    plt.xlabel(r'$\phi$(deg)', fontsize = 20)
    plt.ylabel(r'$m_E(\phi)$', fontsize = 20)
    plt.grid()
    plt.legend(loc = 0, frameon = False, numpoints = 1)
    plt.show()
    





  
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
    
