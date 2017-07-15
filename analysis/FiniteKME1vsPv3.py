import sys
from scipy.optimize import fsolve
import numpy as np
from scipy.integrate import quad
from scipy.special import erfc, erfcinv
import pylab as plt
import sys
basefolder = "/homecentral/srao/Documents/code/mypybox"
sys.path.append(basefolder)
sys.path.append(basefolder + "/utils")
from Print2Pdf import Print2Pdf

def alphaE(m0, meanFieldRates):
    return cFF * JE0**2 * m0 + JEE**2 * meanFieldRates[0] + JEI**2 * meanFieldRates[1]

def MyErfc(z):
    # normal CDF
    return 0.5 * erfc(z / np.sqrt(2.0))

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

def HPrime(z):
    return -1.0 * np.exp(-0.50 * z**2) / np.sqrt(2.0 * np.pi)

def LoadFrRates(N, K, p, m_ext, trialNo, T, gamma = 0, IF_VERBOSE = False, IF_ONLY_E = False):
    fldrName = "./data/N%sK%sm0%sp%sgamma%sT%s/"%(N, K, int(m_ext * 1e3), int(p*10), int(gamma), int(T * 1e-3))
    filename = '/meanrates_theta0_tr%s.txt'%(trialNo)    
    if IF_VERBOSE:
        print fldrName
    if IF_VERBOSE:
        print filename
    out = np.loadtxt(fldrName + filename)
    if IF_ONLY_E:
        out = out[:N]
    return out

def NonSmoothedMeanM1(N, K, p, m_ext, nTr, T, gamma = 0):
    m1 = []
    winSize = int(float(N) / 100.0)
    for n in range(nTr):
        try:
            m = LoadFrRates(N, K, p, m_ext, n, T)
            m = m[:N]
            m1.append(M1Component2(m))
        except IOError:
            print 'file not found'
    return np.mean(m1), np.std(m1), np.std(m1) / np.sqrt(float(nTr))

def NonSmoothedMeanM0(N, K, p, m_ext, nTr, T, gamma = 0):
    m0 = []
    winSize = int(float(N) / 100.0)
    for n in range(nTr):
        try:
            m = LoadFrRates(N, K, p, m_ext, n, T)
            m = m[:N]
            m0.append(np.mean(m))
        except IOError:
            print 'file not found'
    return np.mean(m0), np.std(m0), np.std(m0) / np.sqrt(float(nTr))

def M1Component2(x):
    dPhi = np.pi / len(x)
    return  2.0 * np.absolute(np.dot(x, np.exp(-2.0j * np.arange(len(x)) * dPhi))) / len(x)

def M1vsp(pList, m_ext, nTr = 10, T = 1000, NList = [10000], KList=[1000]):
    m1 = []
    plt.figure()
    mPrediction = np.load('./data/analysisData/p_vs_mI1_m0%s.npy'%(int(m_ext*1e3)))
    plt.plot(mPrediction[0, :], mPrediction[1, :] , color = 'k', linestyle = '--', lw = 1, label = 'E, MF')
    IF_SMOOTHED = False
    for k, kK in enumerate(KList):
        for n, nNI in enumerate(NList):
            m1 = []
            semM1 = []
            if(n == 0):
                print '-'*25
                print "N = ", nNI, " K = ", kK
            for iip, p in enumerate(pList):
                tmpM1, dummy, tmpsem = NonSmoothedMeanM1(nNI, kK, p, m_ext, nTr, T)
                m1.append(tmpM1)
                semM1.append(tmpsem)
            m1 = np.array(m1)
            validPIdx = ~np.isnan(m1)
            semM1 = np.array(semM1)
            print m1[validPIdx].shape, np.sum(validPIdx)
            if(m1.size > 0):
                (_, caps, _) = plt.errorbar(pList[validPIdx], m1[validPIdx], fmt = '.-', yerr = semM1[validPIdx], lw = 1.0, elinewidth=1, markersize = 1, label = r'$N = %s, K = %s$'%(nNI, kK))
            for cap in caps:
                cap.set_markeredgewidth(0.2)
            print '-'*25            
    plt.xlabel(r'$p$')    
    plt.ylabel(r'$\left[ m^{(1)} \right]_J$') # + ' avg over 10 realizations')
    plt.title(r"$m_0 = %s, nonsmoothed$"%(m0))
    # plt.legend(loc = 0, frameon = False, numpoints = 1)

def M0vsp(pList, m_ext, nTr = 10, T = 1000, NList = [10000], KList=[1000]):
    m1 = []
    plt.figure()
    for k, kK in enumerate(KList):
        for n, nNI in enumerate(NList):
            m1 = []
            semM1 = []
            for iip, p in enumerate(pList):
                tmpM1, dummy, tmpsem = NonSmoothedMeanM0(nNI, kK, p, m_ext, nTr, T)
                m1.append(tmpM1)
                semM1.append(tmpsem)
            m1 = np.array(m1)
            validPIdx = ~np.isnan(m1)
            semM1 = np.array(semM1)
            print m1[validPIdx].shape, np.sum(validPIdx)
            (_, caps, _) = plt.errorbar(pList[validPIdx], m1[validPIdx], fmt = '.-', yerr = semM1[validPIdx], lw = 1.0, elinewidth=1, markersize = 1, label = r'$N = %s, K = %s$'%(nNI, kK))
            for cap in caps:
                cap.set_markeredgewidth(0.2)
    plt.xlabel(r'$p$')    
    plt.ylabel(r'$\left[ m^{(1)} \right]_J$') # + ' avg over 10 realizations')
    plt.title(r"$m_0 = %s, nonsmoothed$"%(int(m0*1e3)))
    plt.legend(loc = 0, frameon = False, numpoints = 1)

def Pcritical(uE0, alpha, JEE):
    b = (1 - uE0) / np.sqrt(alpha)
    C0 = np.exp(-b**2 / 2.0) / np.sqrt(2.0 * np.pi)
    num = np.sqrt(alpha)
    denom = C0 * JEE
    denum2 = C0 * (JEE + JEE**2 / (2.0 * alpha * np.sqrt(K)))
    print "new corrected p_c = ", num / denum2
    return num / denom

def QFuncInv(z):
    return np.sqrt(2.0) * erfcinv(z * 2.0)
    
def PcriticalFiniteK(uE0, alphaE0, JEE, JEI, hPrime, K, A):
    b = (1 - uE0) / np.sqrt(alphaE0)
    # C0 = np.exp(-b**2 / 2.0) / np.sqrt(2.0 * np.pi)
    num = -1.0 * np.sqrt(alphaE0)
    denomMatlab = hPrime * (JEE + JEE**2 / (2.0 * alphaE0 * np.sqrt(float(K))))
    denom = hPrime * (JEE + JEE**2 * A / (2.0  * np.sqrt(float(K) * alphaE0)))
    print "New corrected p_c = ", num / denom
    print 'denom matlab = ' , denomMatlab

    # alphaKInf = JEE**2 * meanFieldRates[0] + JEI**2 * meanFieldRates[1]
    # num2 = -1.0 * np.sqrt(alphaKInf)
    # hPrimeKInf = HPrime(QFuncInv(meanFieldRates[0]))
    # denom2 = hPrimeKInf * (JEE + JEE**2 * QFuncInv(meanFieldRates[0]) / (2.0  * np.sqrt(float(K) * alphaKInf)))
    # print "p_cKinf = ", num2/denom2
    # raise SystemExit
    return num / denom

def UE0GivenUE1(u0Guess, *args):
    # solve for u_E(0) given U_E(1)
    u0 = u0Guess
    jee, alpha, u1, meanFieldRate = args
    funcME0 = lambda phi: (1.0 / np.pi) * MyErfc((1.0 - u0 - u1 * np.cos(2.0 * phi)) / np.sqrt(alpha))      return meanFieldRate - quad(funcME0, 0, np.pi)[0]

def Eq0(u0Guess, *args):
    # solve for u_E(0)
    u0 = u0Guess
    jee, alpha, me0 = args
    return me0 - MyErfc((1.0 - u0) / np.sqrt(alpha))

def EqFiniteK(u0Guess, *args):
    # solve for u_E(0), u_I(0), u_E(0) at p
    uE0, uI0, mE1 = u0Guess
    K, p, m_ext = args
    detJab = np.linalg.det(Jab)
    sqrtK = np.sqrt(K)
    uTildeE0 = uE0 / sqrtK - JE0 * m_ext * cFF
    uTildeI0 = uI0 / sqrtK - JI0 * m_ext * cFF
    mE0 = ( JII * uTildeE0 - JEI * uTildeI0) / detJab
    mI0 = (-JIE * uTildeE0 + JEE * uTildeI0) / detJab
    uE1 = p * JEE * mE1    
    aE0 = cFF * JE0**2 * m_ext + JEE**2 * mE0 + JEI**2 * mI0
    aI0 = cFF * JI0**2 * m_ext + JIE**2 * mE0 + JII**2 * mI0
    aE1 = p * JEE**2 * mE1 / sqrtK
    aI1 = p * JEI**2 * mE1 / sqrtK    
    if(aE1 > aE0):
        aE1 = aE0 * 0.1
    AEofPhi = lambda phi: aE0 + aE1 * np.cos(2.0 * phi)
    AIofPhi = lambda phi: aI0 + aI1 * np.cos(2.0 * phi)    
    funcME0 = lambda phi: (1.0 / np.pi) * MyErfc((1.0 - uE0 - uE1 * np.cos(2.0 * phi)) / np.sqrt(AEofPhi(phi)))
    funcMI0 = lambda phi: (1.0 / np.pi) * MyErfc((1.0 - uI0) / np.sqrt(AIofPhi(phi)))
    funcME1 = lambda phi: (2.0 * np.cos(2.0 * phi) / np.pi) * MyErfc((1.0 - uE0 - uE1 * np.cos(2.0 * phi)) / np.sqrt(AEofPhi(phi)))
    # print uE0, uI0, AEofPhi(0)
    a = mE0 - quad(funcME0, 0, np.pi, limit = 200, epsabs = 1e-10)[0]
    b = mI0 - np.abs(quad(funcMI0, 0, np.pi, limit = 200, epsabs = 1e-10)[0])
    # b = mI0 - (1.0 / np.pi) * MyErfc((1.0 - uI0) / np.sqrt(aI0))
    c = mE1 - np.abs(quad(funcME1, 0, np.pi, limit = 200, epsabs = 1e-10)[0])
    # print mE0, mE1, mI0
    return (a, b, c)

def EqFiniteKFunc(u0Guess, args, Jab, JE0, JI0):
    # solve for u_E(0), u_I(0), u_E(0) at p
    uE0, uI0, mE1 = u0Guess
    K, p, m_ext = args
    detJab = np.linalg.det(Jab)
    sqrtK = np.sqrt(K)
    uTildeE0 = uE0 / sqrtK - JE0 * m_ext
    uTildeI0 = uI0 / sqrtK - JI0 * m_ext
    JEE = Jab[0][0]
    JEI = Jab[0][1]
    JIE = Jab[1][0]
    JII = Jab[1][1]
    mE0 = ( JII * uTildeE0 - JEI * uTildeI0) / detJab
    mI0 = (-JIE * uTildeE0 + JEE * uTildeI0) / detJab
    uE1 = p * JEE * mE1    
    aE0 = JEE**2 * mE0 + JEI**2 * mI0
    aI0 = JIE**2 * mE0 + JII**2 * mI0
    aE1 = p * JEE**2 * mE1 / sqrtK
    if(aE1 > aE0):
        aE1 = aE0 * 0.1
    AEofPhi = lambda phi: aE0 + aE1 * np.cos(2.0 * phi)
    funcME0 = lambda phi: (1.0 / np.pi) * MyErfc((1.0 - uE0 - uE1 * np.cos(2.0 * phi)) / np.sqrt(AEofPhi(phi)))
    funcMI0 = lambda phi: (1.0 / np.pi) * MyErfc((1.0 - uI0) / np.sqrt(aI0))
    funcME1 = lambda phi: (2.0 * np.cos(2.0 * phi) / np.pi) * MyErfc((1.0 - uE0 - uE1 * np.cos(2.0 * phi)) / np.sqrt(AEofPhi(phi)))
    a = quad(funcME0, 0, np.pi, limit = 200, epsabs = 1e-18)[0]
    b = np.abs(quad(funcMI0, 0, np.pi, limit = 200, epsabs = 1e-18)[0])    
    # b = (1.0 / np.pi) * MyErfc((1.0 - uI0) / np.sqrt(aI0))
    c = np.abs(quad(funcME1, 0, np.pi, limit = 200, epsabs = 1e-18)[0])
    uE0 = sqrtK * (JEE * a + JEI * b + JE0 * m_ext)
    uI0 = sqrtK * (JIE * a + JII * b + JI0 * m_ext)
    # return (uE0, uI0, c)
    return (a, b, c)

def IterativeSolution(nIter, u0Guess, args, absTol = 1e-24, IF_VERBOSE = False):
    uE0, uI0, mE1 = u0Guess
    K, p, m_ext = args
    iterCnt = 0
    old0 = uE0
    old1 = uI0
    old2 = mE1
    tolerance = absTol
    detJab = np.linalg.det(Jab)
    sqrtK = np.sqrt(K)
    if IF_VERBOSE:
        print '-'*26        
        print 'initial values =', uE0, uI0, mE1
        print 'abs diff:'
    while(iterCnt < nIter):
        new0, new1, new2 =  EqFiniteKFunc(u0Guess, args)
        diff0 = np.abs(new0 - old0)
        diff1 = np.abs(new1 - old1)
        diff2 = np.abs(new2 - old2)
        old0 = new0; old1 = new1; old2 = new2;
        if IF_VERBOSE:
            print diff0, diff1, diff2
        if(diff0 < tolerance and diff1 < tolerance and diff2 < tolerance):
            if IF_VERBOSE:
                print 'converged, #iter = ', iterCnt
                print 'uE0, uI0, mE1 = ', new0, new1, new2
            uTildeE0 = new0 / sqrtK - JE0 * m0
            uTildeI0 = new1 / sqrtK - JI0 * m0
            mE0 = ( JII * uTildeE0 - JEI * uTildeI0) / detJab
            mI0 = (-JIE * uTildeE0 + JEE * uTildeI0) / detJab
            if IF_VERBOSE:
                print 'mE0, mI0, mE1 = ',  mE0, mE1, new2
            return new0, new1, new2
        iterCnt = iterCnt + 1
    return np.nan, np.nan, np.nan

def SolveForPCritical(m0List):
    nPoints = 10
    uE1List = np.linspace(1e-6, 1.0, nPoints)
    for l, lM0 in enumerate(m0List):
        meanFieldRates = -1.0 * np.dot(np.linalg.inv(Jab), Ea) * lM0
        alpha_E = alphaE(lM0, meanFieldRates)
        #============================================ 
        uE0Solution = fsolve(Eq0, [0.1], args = (JEE, alphaE(lM0, meanFieldRates), meanFieldRates[0]), full_output = 1)    
        PrintError(uE0Solution, 'u_E(0)')
        uE0 = uE0Solution[0]
        pCAnalytic = Pcritical(uE0, alpha_E, JEE)
        PcriticalCorrected(uE0, alpha_E, JEE);
        # raise SystemExit
        #------------------------------------------
        uE0List = []
        mE1List = []
        pList = []
        for i, iUE1 in enumerate(uE1List):
            uE0Sol = fsolve(UE0GivenUE1, [0.1], args = (JEE, alpha_E, iUE1, meanFieldRates[0]), full_output = True, xtol = 1e-6)
            if(uE0Sol[-2] == 1): # .i.e IF SOLUTION CONVERGED
                FuncME0 = lambda phi: (1.0 / np.pi) * MyErfc((1.0 - uE0Sol[0][0] - iUE1 * np.cos(2.0 * phi)) / np.sqrt(alpha_E))
                # PrintError(uE0Sol, 'U_E(0)[U_E(1)]')            
                print uE0Sol[0][0], iUE1, quad(FuncME0, 0, np.pi)[0], uE0Sol[-1], uE0Sol[-2]
                uE0List.append((uE0Sol[0][0], iUE1))
                FuncME1 = lambda phi: (2.0 / np.pi) * np.cos(2.0 * phi) * MyErfc((1.0 - uE0Sol[0][0] - iUE1 * np.cos(2.0 * phi)) / np.sqrt(alpha_E))
                me1 = quad(FuncME1, 0, np.pi)[0]
                if(me1 <= meanFieldRates[0] and iUE1 > 0):
                    pList.append(iUE1 / (JEE * me1))
                    mE1List.append(me1)
        if(lM0 == m0List[-1]):
            plt.plot(pList, mE1List, 'bo-', markersize = 4, label = 'numeric')
            plt.vlines(pCAnalytic, 0, meanFieldRates[0], 'k', label = 'analytic')
            plt.legend(loc = 0, frameon = False, numpoints = 1)
            plt.xlabel('p')
            plt.ylabel(r'$m_E^{(1)}$')
            # plt.xlim(3.5, 5.5)
            # plt.ylim(0, 0.3)
            plt.text(pList[-1], mE1List[-1], '%s'%(lM0))            
#            plt.show()
        else:
            plt.plot(pList, mE1List, 'bo-', markersize = 4)
            plt.vlines(pCAnalytic, 0, meanFieldRates[0], 'k')
        if(l == 0):
            plt.text(pList[-1], mE1List[-1], r'$m_0 = %s$'%(lM0))
        else:
            plt.text(pList[-1], mE1List[-1], '%s'%(lM0))
    print "--" * 25
    print "m0 = ", m0List[0], "p_critical = ", pCAnalytic[0]
    print "--" * 24    
    np.save('./data/p_vs_mI1_m0%s'%(int(lM0 * 1e3)), [pList, mE1List, np.ones((len(pList), )) * meanFieldRates[0], np.ones((len(pList), )) * meanFieldRates[1]])            

# def CirConvolve(signal, windowSize):
#     ker = np.concatenate((np.ones((windowSize, )), np.zeros((signal.size - windowSize, ))))
#     return np.real(np.fft.ifft( np.fft.fft(signal)*np.fft.fft(ker))) * (1.0 / float(windowSize))
def CirConvolve(x, winSize):
    print x.shape, winSize
    window = np.ones((winSize, ))
    return convolve(x, winSize, mode = 'wrap')

def PolyFit(pVal, me1, pCAnalytic):
    x = pVal[pVal > pCAnalytic + 0.1]
    y = me1[pVal > pCAnalytic + 0.1]
    # polyCoeff = np.polyfit(np.sqrt(x)

def integrand(x):
    uE0, uI0, mE1 = u0Guess
    K, p, m_ext = args

    uE0 = x[0]
    uI0 = x[1]
    mE1 = x[2]
    beta  = x[3]
    gamma = x[4]
    phi   = x[5]
    k = 1.
    T = 1.
    ww = w(r, theta, phi, alpha, beta, gamma)
    return (math.exp(-ww/(k*T)) - 1.)*r*r*math.sin(beta)*math.sin(theta)
    
def monte_carlo_two(function, a, b, numberOfseconds):
      time_0 = time.clock()
      sumOfpoints = 0
      numberOftrials = 0
      while abs(time_0 - time.clock()) <= numberOfseconds:
          x = Decimal(random.uniform(a, b))
          sumOfpoints += Decimal(function(x))
          numberOftrials += 1
          if numberOftrials == 0:
              return 0
          average = sumOfpoints / Decimal(numberOftrials)
      return average * Decimal(b - a)

if __name__ == "__main__":
    JEE = 1.0
    JIE = 1.0
    JEI = -1.5 
    JII = -1.0
    cFF = 0.2
    JE0 = 2.0 / cFF
    JI0 = 1.0 / cFF
    Jab = np.array([[JEE, JEI],
                    [JIE, JII]])
    Ea = np.array([JE0, JI0])

    CheckBalConditions(JEE, JEI, JIE, JII, JE0, JI0)
    gamma = 0.0
    m0 = 0.075
    p = 7
    K = int(sys.argv[2])
    simTime = int(sys.argv[3])
    solverType = 'fsolve'
    if(len(sys.argv) > 3):
        solverTypeStr = sys.argv[4]
        if(solverTypeStr == 'fpi'):
           solverType = 'fixed-point-iter'
    meanFieldRates = -1.0 * np.dot(np.linalg.inv(Jab), Ea) * m0
    print "mf rates == ", meanFieldRates
    print "alphaE mf == ", JEE**2 * meanFieldRates[0] + JEI**2 * meanFieldRates[1]
    
    aE0 = JEE**2 * meanFieldRates[0]+ JEI**2 * meanFieldRates[1]
    aI0 = JIE**2 * meanFieldRates[0]+ JII**2 * meanFieldRates[1]        
    sqrtK = np.sqrt(K)
    uE0Sol = fsolve(Eq0, [1e-2], args = (JEE, aE0, meanFieldRates[0]), full_output = True, xtol = 1e-6)
    uI0Sol = fsolve(Eq0, [1e-2], args = (JII, aI0, meanFieldRates[1]), full_output = True, xtol = 1e-6)
    u0Guess = [uE0Sol[0][0], uI0Sol[0][0], 0.2613249636852059 ]

    # u0Guess = [uE0Sol[0][0], uI0Sol[0][0], 1e-6]    
    print "initial guess = ", u0Guess

    uE0Sol = fsolve(EqFiniteK, u0Guess, args = (K, 1, m0), full_output = True, xtol = 1e-12, epsfcn = 1e-12)
    PrintError(uE0Sol, 'u_E(0),  u_I(0), m_E(1)]')
    # print "iterating...  "
    # nIterations = 2000
    # IterativeSolution(nIterations, uE0Sol[0], (K, p, m0), IF_VERBOSE = True)
    # raise SystemExit
    uE0, uI0, mE1 = uE0Sol[0]
    detJab = np.linalg.det(Jab)
    sqrtK = np.sqrt(K)
    uTildeE0 = uE0 / sqrtK - JE0 * m0 * cFF
    uTildeI0 = uI0 / sqrtK - JI0 * m0 * cFF
    mE0 = ( JII * uTildeE0 - JEI * uTildeI0) / detJab
    mI0 = (-JIE * uTildeE0 + JEE * uTildeI0) / detJab
    aE0 = cFF * JE0**2 * m0 + JEE**2 * mE0 + JEI**2 * mI0
    aI0 = cFF * JI0**2 * m0 + JIE**2 * mE0 + JII**2 * mI0
    print "-"*26
    hPrimeIn =  HPrime((1.0 - uE0)/ np.sqrt(aE0))
    pCriticalAtK = PcriticalFiniteK(uE0, aE0, JEE, JEI, hPrimeIn, K, (1.0 - uE0)/ np.sqrt(aE0))
    print "m0 = ", m0
    print 'alphaE0 = ', aE0, 'hprime = ', hPrimeIn
    # raise SystemExit
    p = float(sys.argv[1])

    print 'x'*50
    print "rates at K = ", K, "me0, mi0, me1"
    print mE0, mI0, mE1

    raise SystemExit
    np.save('./data/pCritical_m0%s_K%s'%(int(m0 * 1e3), K), pCriticalAtK)

    
    pListSim = np.array([0, .5, 1, 2, 2.5, 2.6, 2.7, 2.8]) #, 2.9, 3.0])
    pListSim = np.arange(0, 4.5, 0.1)
    pStart = pCriticalAtK + 2
    pList = np.linspace(pStart, pCriticalAtK - 1, 5000)
    # pList = np.linspace(pStart, pCriticalAtK - 1, 10000)    
    # pList = np.linspace(10, 4, 8000)    
    # pList = np.linspace(7, 3, 40000)        
    # pList = np.logspace(np.log10(pCritticalatK - 0.2), np.log10(pCritticalatK + 2), 5000)
    nPoints = len(pList)
    print 'nPoints = ', nPoints
    mE1SolList = []
    mE0SolList = []
    mI0SolList = []        
    pValidList = []
    uE0SolList = []
    uI0SolList = []    
    # simRates = np.loadtxt('
    # M1vsp(pListSim, m0, nTr = 10, T = simTime, NList = [10000, 20000, 40000], KList=[K])
    # uE0Sol = fsolve(EqFiniteK, [.943, .5, 0.8 / pStart], args = (K, pStart, m0), full_output = True, xtol = 1e-12, epsfcn = 1e-12)
    uE0Sol = fsolve(EqFiniteK, u0Guess, args = (K, pStart, m0), full_output = True, xtol = 1e-12, epsfcn = 1e-12)    
    PrintError(uE0Sol, 'u_E(0),  u_I(0), m_E(1)] at p = %s'%(pStart))
    # u0Guess = [uE0Sol[0][0], uI0Sol[0][0], 0.1606553588900383]    
    print "numerics ..."    
    for ip, p in enumerate(pList):
        if ip == 0:

            uE0Sol = fsolve(EqFiniteK, u0Guess, args = (K, p, m0), full_output = True, xtol = 1e-10, epsfcn = 1e-10)
            
            # if solverType == 'fixed-point-iter':
            #     uE0Sol = IterativeSolution(nIterations, uE0Sol[0], (K, p, m0))
            # else:
            #     uE0Sol = fsolve(EqFiniteK, u0Guess, args = (K, p, m0), full_output = True, xtol = 1e-10, epsfcn = 1e-10)
        else:
            if not (len(uE0SolList) == 0 or len(uI0SolList) == 0 or len(mE1SolList) == 0):
               u0Guess = [uE0SolList[-1], uI0SolList[-1], np.abs(mE1SolList[-1]) * 1e-1]
#              u0Guess = [uE0SolList[-1], uI0SolList[-1], 0.1606553588900383]
            uE0Sol = fsolve(EqFiniteK, uE0Sol[0], args = (K, p, m0), full_output = True, xtol = 1e-10, epsfcn = 1e-10)
            # if solverType == 'fixed-point-iter':
            #     uE0Sol = IterativeSolution(nIterations, [uE0SolList[-1], uI0SolList[-1], np.abs(mE1SolList[-1]) * 1e-1], (K, p, m0))
            # else:
            #     uE0Sol = fsolve(EqFiniteK, [uE0SolList[-1], uI0SolList[-1], np.abs(mE1SolList[-1]) * 1e-1], args = (K, p, m0), full_output = True, xtol = 1e-10, epsfcn = 1e-10)

        if solverType == 'fixed-point-iter':
            uE0Sol = IterativeSolution(nIterations, uE0Sol[0], (K, p, m0), IF_VERBOSE = True)
            uE0, uI0, mE1 = uE0Sol
            IF_CONVERGED = ~np.isnan(uE0) and ~np.isnan(uI0) and ~np.isnan(mE1)
        else:
            uE0, uI0, mE1 = uE0Sol[0]
            IF_CONVERGED = uE0Sol[-2] == 1
        uTildeE0 = uE0 / sqrtK - JE0 * m0
        uTildeI0 = uI0 / sqrtK - JI0 * m0
        mE0 = ( JII * uTildeE0 - JEI * uTildeI0) / detJab
        mI0 = (-JIE * uTildeE0 + JEE * uTildeI0) / detJab
        aE0 = JEE**2 * mE0 + JEI**2 * mI0
        aI0 = JIE**2 * mE0 + JII**2 * mI0    
        if(IF_CONVERGED): # .i.e IF SOLUTION CONVERGED
            uE0SolList.append(uE0)
            uI0SolList.append(uI0)            
            mE0SolList.append(mE0)
            mI0SolList.append(mI0)
            mE1SolList.append(mE1)            
            pValidList.append(p) 
    print "done"   
    plt.plot(pValidList, np.array(mE1SolList, dtype = 'float'), 'c.', label = 'corrected', lw = 1.5)
    np.save('./data/mE1_m0%s_K%s'%(int(m0 * 1e3), K), [pValidList, np.array(mE1SolList, dtype = 'float')])

#    np.save('./data/analysisData/p_vs_mE1_m0%s'%(int(m0*1e3)), [pValidList, np.array(mE1SolList, dtype = 'float'), np.ones((len(pValidList, )) * , np.ones((len(pofME), )) * )    
    #np.save('./data/analysisData/mE1_m0%s_K%s'%(int(m0 * 1e3), K), [pValidList, np.array(mE1SolList, dtype = 'float')])    

    plt.figure()
    plt.plot(pValidList, np.array(mE1SolList, dtype = 'float')**2, 'c.', label = 'corrected', lw = 1.5)
    plt.xlabel(r'$p$')
    plt.ylabel(r'$\left( m_E^{(1)} \right)^2$')


    ymax = plt.ylim()
    plt.vlines(pCriticalAtK, 0, ymax[1], 'k')

    # raise SystemExit
    
    print len(pValidList)
    if(len(pValidList) > 0):
        try:
            plt.legend(loc = 0, frameon = False, numpoints = 1)
        except IndexError:
            print 'no smimulation points found'
    plt.title(r"$m_0 = %s$"%(m0))    
    figname = 'mE1_m0%s_K%s.png'%(int(m0 * 1e3), K)
    print 'saving as', figname
    plt.savefig('./figs/' + figname)# m1_vs_N_m0%s_p%s.png'%(int(m_ext * 1e3), int(p*100)))
    plt.ylabel(r'$\left[ m^{(0)} \right]_J$') # + ' avg over 10 realizations')
    plt.figure()
    plt.plot(pValidList, np.array(mE0SolList, dtype = 'float'), 'k-', label = 'me0', lw = 1.5)
    plt.plot(pValidList, np.array(mI0SolList, dtype = 'float'), 'r-', label = 'mi0', lw = 1.5)
    plt.legend()
    
    plt.show()
    raise SystemExit

    M0vsp(pListSim, m0, nTr = 10, T = 1000, NList = [10000], KList=[1000])
    mPrediction = np.load('./data/analysisData/p_vs_mI1_m0%s.npy'%(int(m0*1e3)))
    plt.hlines(meanFieldRates[0], 0, pList[-1], 'k') #,  plot(mPrediction[0, :], mPrediction[2, :] , color = 'k', linestyle = '--', lw = 1, label = 'E, MF')
    plt.plot(pValidList, mE0SolList, 'go-', label = 'corrected')    
    plt.xlabel(r'$p$')    
    plt.ylabel(r'$\left[ m^{(0)} \right]_J$') # + ' avg over 10 realizations')
    plt.title(r"$m_0 = %s, nonsmoothed$"%(m0))
    plt.legend(loc = 0, frameon = False, numpoints = 1)
    figname = 'nonsmoothed_mE0_m0%s.png'%(int(m0 * 1e3))
    print 'saving as', figname
    plt.savefig('./figs/' + figname)# m1_vs_N_m0%s_p%s.png'%(int(m_ext * 1e3), int(p*100)))
    
    plt.show()
    
    raise SystemExit
    m0List = [0.075, 0.1, 0.15, 0.175, .2, .3]
    m0List = [0.075]
    SolveForPCritical(m0List)
 
    figFolder = '' 
    figname = './p_vs_mE1' # '/p%s'%(p)
    paperSize = [5, 4.0]
    figFormat = 'png'
    axPosition = [0.17, 0.15, .8, 0.8]
    print figFolder, figname
    Print2Pdf(plt.gcf(),  figFolder + figname,  paperSize, figFormat=figFormat, labelFontsize = 12, tickFontsize=10, titleSize = 10.0, IF_ADJUST_POSITION = True, axPosition = axPosition)
    plt.show()



   
