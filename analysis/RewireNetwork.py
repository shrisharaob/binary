''' USAGE: [dbName, NetworkType, K, NE, NI, thetaSig] = DefaultArgs(sys.argv[1:], ['', 'ori', 1000, 10000, 10000, 0.5]) '''

basefolder = "/homecentral/srao/Documents/code/mypybox"
import numpy as np
import code, sys, os, ipdb
import pylab as plt
sys.path.append(basefolder)
import Keyboard as kb
from multiprocessing import Pool
from functools import partial 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
sys.path.append(basefolder + "/nda/spkStats")
sys.path.append(basefolder + "/utils")
from DefaultArgs import DefaultArgs
from reportfig import ReportFig
from Print2Pdf import Print2Pdf
import GetPO


# probFunc = lambda thetaDiff, sig: 1.0 / (sig * np.sqrt(2.0*np.pi)) * np.exp(-(np.sin((np.pi/180.0)*(thetaDiff)))**2 / (2.0*sig**2))
# out = probFunc(thetaDiff, sig)
# out[thetaDiff >= np.pi / 2.0] = (1.0/(sig * np.sqrt(2.0 * np.pi))) - out[thetaDiff >= np.pi/2.0] 
# return 1.0 + np.cos(2.0 * thetaDiff)
probFunc = lambda thetaDiff, sig: 1.0 / (sig * np.sqrt(2.0*np.pi)) * np.exp(-(np.sin((thetaDiff)))**2 / (2.0*sig**2))
def ConProbFunc(thetaDiff, sig):
    #thetaDiff in radians
    return probFunc(thetaDiff, sig)

# os.system("gcc -fPIC -o gensparsevec.so -shared GenSparseVec.c") 
# #os.system("gcc -g -ggdb -o gensparsevec.so -shared GenSparseVec.c")
# mycfunc = np.ctypeslib.load_library('gensparsevec.so', '.') # use gcc -o gensparsevec.so -shared GenSparseVec.c
# mycfunc.GenSparseMat.restype = None
# mycfunc.GenSparseMat.argtypes = [np.ctypeslib.ndpointer(np.int32, flags = 'aligned, contiguous'),
#                                  np.ctypeslib.c_intp,
#                                  np.ctypeslib.c_intp,
#                                  np.ctypeslib.ndpointer(np.int32, flags = 'aligned, contiguous, writeable'),
#                                  np.ctypeslib.ndpointer(np.int32, flags = 'aligned, contiguous, writeable'),
#                                  np.ctypeslib.ndpointer(np.int32, flags = 'aligned, contiguous, writeable')]

def GenSparseMat(convec, rows, clmns, sparseVec, idxvec, nPostNeurons):
    requires = ['CONTIGUOUS', 'ALIGNED']
    convec = np.require(convec, np.int32, requires)
    sparseVec = np.require(sparseVec, np.int32, requires)
    idxvec = np.require(idxvec, np.int32, requires)
    nPostNeurons = np.require(nPostNeurons, np.int32, requires);
    mycfunc.GenSparseMat(convec, rows, clmns, sparseVec, idxvec, nPostNeurons)
    fpsparsevec = open('sparseConVec.dat', 'wb')
    sparseVec.tofile(fpsparsevec)
    fpsparsevec.close()
    fpIdxVec = open('idxVec.dat', 'wb')
    idxvec.tofile(fpIdxVec)
    fpIdxVec.close()
    fpNpostNeurons = open('nPostNeurons.dat', 'wb')
    nPostNeurons.tofile(fpNpostNeurons)
    fpNpostNeurons.close()

def GenerateFixedInDegreeMat(cprob, K, NE, NI):
    N = NE + NI
    print 'generating fixed K mat'
    for j in range(N):
        preNeurons = np.random.choice(NE, K, replace = False) # 2 for two populations
        cprob[preNeurons, j] = 1.0 # E --> E
        preNeurons = np.random.choice(np.arange(NE, N), K, replace = False) 
        cprob[preNeurons, j] = 1.0 # I -> E
        # cprob[i, preNeurons] = 1.0
    print cprob
    return cprob
            
def Convert2InDegree(idxVec, nPostNeurons, sparseConVec):
    N = len(nPostNeurons)
    nPreNeurons = np.zeros((N, ))
    preNeurons = [[] for x in xrange(N)]
    # preIdxVec = np.zeros((N, ))
    for i in range(N):
	iPostNeurons = sparseConVec[idxVec[i] : idxVec[i] + nPostNeurons[i]]
	for j, iPost in enumerate(iPostNeurons):
	    preNeurons[iPost].append(i)
    return preNeurons

def Convert2OutDegree(preNeurons, nPostNeurons):
    # sparseConVec = np.zeros((nPostNeurons.sum(), ))
    N = len(nPostNeurons)    
    sparseConVec = [[] for x in xrange(N)]
    count = 0
    for i in range(N):
	iPreNeurons = preNeurons[i]
	# ipdb.set_trace()
	for kk in iPreNeurons:
	    sparseConVec[kk].append(i)
    return np.asarray(np.hstack(sparseConVec), dtype = np.int32)

def RewireSqrtK(preNeurons, recModulation, po, K, NE):
    N = NE #len(preNeurons)
    recModulation = 1
    print 'nPre befor = ', len(np.hstack(preNeurons))
    for i in range(N):
	iPreNeurons = preNeurons[i]	
	iPreNeurons = np.array(iPreNeurons)
	preLenghtOld = iPreNeurons.size
	iPreNeuronsI = iPreNeurons[iPreNeurons >= N]
	iPreNeurons = iPreNeurons[iPreNeurons < N]
	iK = len(iPreNeurons)
	nLinks2Break = int(np.sqrt(float(iK)))
    	links2Break = np.sort(np.random.choice(range(iK), nLinks2Break, replace = False))
	maxIters = nLinks2Break * 50
	if maxIters > N:
	    maxIters = N
	randomPreNeuron = np.random.permutation(xrange(N)) 
	vidx = np.empty((randomPreNeuron.size, ))
	vidx[:] = True
	for mm in links2Break:
	    vidx = np.logical_and(randomPreNeuron != iPreNeurons[mm], vidx)
	randomPreNeuron = randomPreNeuron[vidx]
	maxIters = randomPreNeuron.size
        iterCount = 0
        nUnique = np.unique(iPreNeurons)
	for j in links2Break:
	    nRewiredCount = 0	
	    while nRewiredCount < 1 and iterCount < maxIters:
		prob = 0.5 * (1 + np.cos(2.0 * (po[i] - po[randomPreNeuron[iterCount]])))
		if(prob >= np.random.rand()):
		    if j > 0:
			if not np.any(iPreNeurons == randomPreNeuron[iterCount]):
			    iPreNeurons[j] = randomPreNeuron[iterCount]
			    nRewiredCount += 1
		iterCount += 1
	preNeurons[i] = iPreNeurons.tolist()
	preNeurons[i].extend(iPreNeuronsI.tolist())
	preLenghtNew = len(preNeurons[i])
	if(preLenghtOld != preLenghtNew):
	    print 'dose not work for neuron#', i
    print 'nPre after = ', len(np.hstack(preNeurons))    
    return preNeurons


def LoadFr(p, gamma, phi, mExt, mExtOne, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False):
    baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
    if nPop == 1:
    	baseFldr = baseFldr + 'onepop/data/rewire/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if nPop == 2:
	if gamma >= .1 or gamma == 0:
	    baseFldr = baseFldr + 'twopop/data/rewire/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
	else:
	    baseFldr = baseFldr + 'twopop/data/rewire/N%sK%s/m0%s/mExtOne%s/p%sgamma/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(T*1e-3), trNo)	    
    if IF_VERBOSE:
    	print baseFldr
    filename = 'meanrates_theta%.6f_tr%s.txt'%(phi, trNo)
    print filename
    return np.loadtxt(baseFldr + filename)


def GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000):
    NE = N
    NI = N
    tc = np.zeros((NE + NI, nPhis))
    tc[:] = np.nan
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    for i, iPhi in enumerate(phis):
	print i, iPhi
	try:
	    if i == 0:
		print 'loading from fldr: ',	    
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, trNo, T, NE, K, nPop, IF_VERBOSE = True)
	    else:
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, trNo, T, NE, K, nPop, IF_VERBOSE = False)
	    if(len(fr) == 1):
		if(np.isnan(fr)):
		    print 'file not found!'
	    tc[:, i] = fr
	except IOError:
	    print 'file not found!'
	    raise SystemExit
    return tc


def GetPhase(firingRate, atTheta, IF_IN_RANGE = False):
    out = np.nan
    zk = np.dot(firingRate, np.exp(2j * atTheta * np.pi / 180))
    out = np.angle(zk) * 180.0 / np.pi
    if IF_IN_RANGE:
	if(out < 0):
	    out += 360
    return out * 0.5

def POofPopulation(tc, theta = np.arange(0.0, 180.0, 22.5), IF_IN_RANGE = False):
    # return value in degrees
    nNeurons, nAngles = tc.shape
    theta = np.linspace(0, 180, nAngles, endpoint = False)
    po = np.zeros((nNeurons, ))
    for kNeuron in np.arange(nNeurons):
        po[kNeuron] = GetPhase(tc[kNeuron, :], theta, IF_IN_RANGE)
    return po 

def GetPOofPop(p, gamma, mExt, mExtOne, nPhis = 8, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_IN_RANGE = True):
    nNeurons = N
    thetas= np.linspace(0, np.pi, nNeurons, endpoint = False)
    tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, trNo, N, K, nPop, T)
    prefferedOri = POofPopulation(tc[:N], IF_IN_RANGE = True) * np.pi / 180.0
    return prefferedOri
    
if __name__ == '__main__':
    # NetworkType : {'uni', 'ori'}, 'uni' is for standard random network, 'ori' is to rewire depending on the distance in ori space
    [trNo, p, gamma, mExt, mExtOne, nPhis, K, NE, NI, thetaSig, thetaSigI, nPop, T] = DefaultArgs(sys.argv[1:], [0, 0, 0, .075, .075, 8, 1000, 10000, 10000, .75, 0.75, 2, 1000])
    NE = int(NE)
    NI = int(NI)
    N = NE
    K = int(K)
    p = float(p)
    trNo = int(trNo)
    mExt = float(mExt)
    mExtOne = float(mExtOne)
    gamma = float(gamma)
    nPop = int(nPop)
    nPhis = int(nPhis)
    thetaSig = float(thetaSig)
    thetaSigI = float(thetaSigI)
    # cprob = np.zeros((NE + NI, NE + NI))
    rootFolder = ''
    baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
    if nPop == 1:
    	baseFldr = baseFldr + 'onepop/data/rewire/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if nPop == 2:
    	if gamma >= .1 or gamma == 0:
    	    baseFldr = baseFldr + 'twopop/data/rewire/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    	if gamma == 0.05:
    	    baseFldr = baseFldr + 'twopop/data/rewire/N%sK%s/m0%s/mExtOne%s/p%sgamma/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(T*1e-3), trNo)	    


    # baseFldr = './'
    # NE = 500
    # NI = 500
    # K = 50

    requires = ['CONTIGUOUS', 'ALIGNED']
    # convec = np.require(convec, np.int32, requires)
    nPostNeurons = np.zeros((NE + NI, ))
    idxvec = np.zeros((NE + NI, ))
    sparseVec = np.zeros((nPostNeurons.sum(), ))
    
    sparseVec = np.require(sparseVec, np.int32, requires)
    idxvec = np.require(idxvec, np.int32, requires)
    nPostNeurons = np.require(nPostNeurons, np.int32, requires);

    # READ
    fpsparsevec = open(baseFldr + '/' + 'sparseConVec.dat', 'rb')
    sparseVec = np.fromfile(fpsparsevec, dtype = np.int32)
    fpsparsevec.close()
    fpIdxVec = open(baseFldr + '/' + 'idxVec.dat', 'rb')
    idxvec = np.fromfile(fpIdxVec, dtype = np.int32)
    fpIdxVec.close()
    fpNpostNeurons = open(baseFldr + '/' + 'nPostNeurons.dat', 'rb')
    nPostNeurons = np.fromfile(fpNpostNeurons, dtype = np.int32)
    fpNpostNeurons.close()    

    nPostE = nPostNeurons[:NE].sum()
    tmp = sparseVec[sparseVec < NE]
    print np.sort(tmp[:10])
    # REWIRE
    po = GetPOofPop(p, gamma, mExt, mExtOne, nPhis, trNo, N, K, nPop, T, IF_IN_RANGE = True)
    preNeurons = Convert2InDegree(idxvec, nPostNeurons, sparseVec)
    # rewiredPreNeurons = RewireSqrtK(preNeurons, p, po, K, NE)
    rewiredPreNeurons = preNeurons
    sparseConVec = Convert2OutDegree(rewiredPreNeurons, nPostNeurons)
    sparseVec = np.require(sparseConVec, np.int32, requires)

    tmp = sparseVec[sparseVec < NE]
    print np.sort(tmp[:10])
    print 'nConnections = ', nPostNeurons.sum()
    print 'sparsevec length = ', sparseVec.size


    trNo += 1 # connectivity written to another folder
    baseFldrNew = rootFolder + '/homecentral/srao/Documents/code/binary/c/'    
    if nPop == 1:
    	baseFldrNew = baseFldrNew + 'onepop/data/rewire/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if nPop == 2:
    	if gamma >= .1 or gamma == 0:
    	    baseFldrNew = baseFldrNew + 'twopop/data/rewire/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    	if gamma == 0.05:
    	    baseFldrNew = baseFldrNew + 'twopop/data/rewire/N%sK%s/m0%s/mExtOne%s/p%sgamma/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(T*1e-3), trNo)	    


    os.system('mkdir -p ' + baseFldrNew)
    os.system('cp ' + baseFldr + '/*.out ' + baseFldrNew)
    os.system('cp ' + baseFldr + '/*FF.dat ' + baseFldrNew)
    

    # WRITE
    print '--'*20
    print 'saving files to: ', baseFldrNew
    fpsparsevec = open(baseFldrNew + '/' + 'sparseConVec.dat', 'wb')
    sparseVec.tofile(fpsparsevec)
    fpsparsevec.close()
    fpIdxVec = open(baseFldrNew + '/' + 'idxVec.dat', 'wb')
    idxvec.tofile(fpIdxVec)
    fpIdxVec.close()
    fpNpostNeurons = open(baseFldrNew + '/' + 'nPostNeurons.dat', 'wb')
    nPostNeurons.tofile(fpNpostNeurons)
    fpNpostNeurons.close()    

