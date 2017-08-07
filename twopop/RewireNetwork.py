''' USAGE: [dbName, NetworkType, K, NE, NI, thetaSig] = DefaultArgs(sys.argv[1:], ['', 'ori', 1000, 10000, 10000, 0.5]) '''

basefolder = "/homecentral/srao/Documents/code/mypybox"
import numpy as np
import code, sys, os
import ipdb
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

def BuildFullMatrix(idxVec, nPostNeurons, sparseConVec, NE):
    conmat = np.zeros((NE, NE))
    for i in range(NE):
	postEofi = sparseConVec[idxVec[i] : idxVec[i] + nPostNeurons[i]]
	postEofi = postEofi[postEofi < NE]
	for j in postEofi:
	    conmat[i, j] = 1
    return np.asarray(conmat, dtype = np.int32)

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

def RewireProbFunc(po0, po1, funcType = 'cos'):
    out = 0
    # if funcType == 'cos':
    # 	out = 0.5 * (1 + np.cos(2.0 * (po0 - po1)))
    if funcType == 'exp':
    	z = 1.0
    	out = 0.5 * (1 - np.cos(2.0 * (po0 - po1)))
    else:
	out = 0.5 * (1 + np.cos(2.0 * (po0 - po1)))
    return out

def GetNLinks(neuronIdx, po, nLinks2Break, iPreNeurons):
    # returns n links that are far in the orientaion space
    po0 = po[neuronIdx]
    maxIters = iPreNeurons.size
    iterCount = 0
    z = 1.0
    kPreNeurons = iPreNeurons.size
    randomPreNeuron = np.random.permutation(iPreNeurons)
    linksRemoved = []
    nLinksCount = 0
    while nLinksCount < nLinks2Break and iterCount < maxIters:
	prob = 1 - z * np.exp(-z * np.abs(po0 - po[randomPreNeuron[iterCount]]))
	if prob < 0 or prob > 1:
	    print 'link removal prob not in range!!!'
	    raise SystemExit
	if prob >= np.random.rand():
	    linksRemoved.append(randomPreNeuron[iterCount])
	    nLinksCount = nLinksCount + 1
	iterCount += 1
    linksRemovedIdx = [np.where(iPreNeurons == lll)[0][0] for lll in linksRemoved]
    #ipdb.set_trace()
    return np.array(linksRemovedIdx)
	
def RewireSqrtK(preNeurons, recModulation, po, K, NE, rewireType, stepNo):
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
	nLinks2Break = int(np.sqrt(float(iK) / float(stepNo)))
	# ipdb.set_trace()
	if rewireType == 'rand':
	    links2Break = np.sort(np.random.choice(range(iK), nLinks2Break, replace = False))
	if rewireType == 'exp':
	    links2Break = GetNLinks(i, po, nLinks2Break, iPreNeurons)
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
		# prob = 0.5 * (1 + np.cos(2.0 * (po[i] - po[randomPreNeuron[iterCount]])))
                prob = RewireProbFunc(po[i], po[randomPreNeuron[iterCount]], 'cos')		
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

def GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2):
    if rewireType == 'rand':
	tag = ''
    if rewireType == 'exp':
	tag = '1'
    rootFolder = ''
    baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
    if nPop == 1:
    	baseFldr = baseFldr + 'onepop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    if nPop == 2:
	if gamma >= .1 or gamma == 0:
	    baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
	else:
	    baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/p%sgamma/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(T*1e-3), trNo)
    print baseFldr	    
    return baseFldr
		
def LoadFr(p, gamma, phi, mExt, mExtOne, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False, rewireType = 'rand'):
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop)
    if IF_VERBOSE:
    	print baseFldr
    filename = 'meanrates_theta%.6f_tr%s.txt'%(phi, trNo)
    print filename
    return np.loadtxt(baseFldr + filename)

def GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000):
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
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, trNo, T, NE, K, nPop, rewireType = rewireType, IF_VERBOSE = True)
	    else:
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, trNo, T, NE, K, nPop, rewireType = rewireType, IF_VERBOSE = False)
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

def GetPOofPop(p, gamma, mExt, mExtOne, rewireType, nPhis = 8, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_IN_RANGE = True):
    nNeurons = N
    thetas= np.linspace(0, np.pi, nNeurons, endpoint = False)
    tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T)
    prefferedOri = POofPopulation(tc[:N], IF_IN_RANGE = True) * np.pi / 180.0
    return prefferedOri

def LoadSparseMat(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop):
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop)
    NE = N
    NI = N
    requires = ['CONTIGUOUS', 'ALIGNED']
    # # convec = np.require(convec, np.int32, requires)
    # nPostNeurons = np.zeros((NE + NI, ))
    # idxvec = np.zeros((NE + NI, ))
    # idxvec = np.require(idxvec, np.int32, requires)
    # nPostNeurons = np.require(nPostNeurons, np.int32, requires);
    # READ
    fpIdxVec = open(baseFldr + 'idxVec.dat', 'rb')
    idxvec = np.fromfile(fpIdxVec, dtype = np.int32)
    fpIdxVec.close()
    fpNpostNeurons = open(baseFldr + 'nPostNeurons.dat', 'rb')
    nPostNeurons = np.fromfile(fpNpostNeurons, dtype = np.int32)
    fpNpostNeurons.close()
    sparseVec = np.zeros((nPostNeurons.sum(), ))
    sparseVec = np.require(sparseVec, np.int32, requires)
    fpsparsevec = open(baseFldr + 'sparseConVec.dat', 'rb')
    sparseVec = np.fromfile(fpsparsevec, dtype = np.int32)
    fpsparsevec.close()
    # requires = ['CONTIGUOUS', 'ALIGNED']
    # # convec = np.require(convec, np.int32, requires)
    # nPostNeurons = np.zeros((NE + NI, ))
    # idxvec = np.zeros((NE + NI, ))
    # idxvec = np.require(idxvec, np.int32, requires)
    # nPostNeurons = np.require(nPostNeurons, np.int32, requires);
    # # READ
    # fpIdxVec = open(baseFldr + 'idxVec.dat', 'rb')
    # idxvec = np.fromfile(fpIdxVec, dtype = np.int32)
    # fpIdxVec.close()
    # fpNpostNeurons = open(baseFldr + 'nPostNeurons.dat', 'rb')
    # nPostNeurons = np.fromfile(fpNpostNeurons, dtype = np.int32)
    # fpNpostNeurons.close()
    # sparseVec = np.zeros((nPostNeurons.sum(), ))
    # sparseVec = np.require(sparseVec, np.int32, requires)
    # fpsparsevec = open(baseFldr + 'sparseConVec.dat', 'rb')
    # sparseVec = np.fromfile(fpsparsevec, dtype = np.int32)
    # fpsparsevec.close()
    return idxvec, nPostNeurons, sparseVec

def ConnectionDistr(p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N = 10000, K = 1000, nPop = 2, T = 1000):
    po = GetPOofPop(p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N, K, nPop, T, IF_IN_RANGE = True)
    idxvec, nPostNeurons, sparseVec = LoadSparseMat(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop)
    preNeurons = Convert2InDegree(idxvec, nPostNeurons, sparseVec)
    z = np.empty((N, ))
    z[:] = np.nan
    for i in range(N):
	tmp = 0
	Ki = 0
	preNeuronsTmp = np.asarray(preNeurons[i], dtype = np.int)
	preNeuronsTmp = preNeuronsTmp[preNeuronsTmp < N]
	# ipdb.set_trace()
	for k in preNeuronsTmp:
	    poi = po[i]
	    pok = po[k]
	    if (~np.isnan(poi)) and  (~np.isnan(pok)):
		Ki += 1
		tmp += np.exp(2j * (poi - pok))
	z[i] = np.abs(tmp) / float(Ki)
    return z

def CompareConDistr(trList, p = 0, gamma = 0, mExt = 0.075, mExtOne = 0.075, rewireType = 'rand', nPhis = 8, N = 10000, K = 1000, nPop = 2, T = 1000):
    plt.figure()
    plt.ion()
    # out = []
    for trNo in trList:
	z = ConnectionDistr(p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N, K, nPop, T)
	# out.append(z)
	plt.hist(z, 50, normed = True, histtype = 'step', label = 'STP: %s'%(trNo))
	plt.show()
    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 8})
    plt.ylabel('Density')
    plt.xlabel(r"$\langle \exp (\phi^{po}_i - \phi^{po}_j) \rangle_j$")
    plt.savefig('./figs/conDistr_p%sg%s.png'%(int(p * 100), int(gamma * 100)))



def GetPODiff(p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N = 10000, K = 1000, nPop = 2, T = 1000):
    po = GetPOofPop(p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N, K, nPop, T, IF_IN_RANGE = True)
    idxvec, nPostNeurons, sparseVec = LoadSparseMat(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop)
    preNeurons = Convert2InDegree(idxvec, nPostNeurons, sparseVec)
    out = []
    for i in range(N):
	preNeuronsTmp = np.asarray(preNeurons[i], dtype = np.int)
	preNeuronsTmp = preNeuronsTmp[preNeuronsTmp < N]
	for k in preNeuronsTmp:
	    # out.append(np.abs(po[i] - po[k]))
	    out.append(np.abs(po[k]))	    
    return np.asarray(out, dtype = np.float)

def CompareDPOvsDist(trList, p = 0, gamma = 0, mExt = 0.075, mExtOne = 0.075, rewireType = 'rand', nPhis = 8, N = 10000, K = 1000, nPop = 2, T = 1000):
    plt.figure()
    plt.ion()
    for trNo in trList:
	z = GetPODiff(p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N, K, nPop, T)
	plt.hist(z, 50, normed = True, histtype = 'step', label = 'STP: %s'%(trNo))
	plt.show()
    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 8})
    plt.ylabel('Density')
    plt.xlabel(r"$\langle \phi^{po}_i - \phi^{po}_j \rangle$")
    plt.savefig('./figs/DeltaPO_vs_conDistr_p%sg%s.png'%(int(p * 100), int(gamma * 100)))



    
  
if __name__ == '__main__':
    # NetworkType : {'uni', 'ori'}, 'uni' is for standard random network, 'ori' is to rewire depending on the distance in ori space
    [trNo, rewireType, p, gamma, mExt, mExtOne, nPhis, K, NE, NI, nPop, T] = DefaultArgs(sys.argv[1:], [0, 'rand', 0, 0, .075, .075, 8, 1000, 10000, 10000,  2, 1000])
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
    T = int(T)
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop)
    IF_TEST = False
    if IF_TEST:
	baseFldr = './'
	NE = 500
	NI = 500
	K = 50
    requires = ['CONTIGUOUS', 'ALIGNED']
    # convec = np.require(convec, np.int32, requires)
    nPostNeurons = np.zeros((NE + NI, ))
    idxvec = np.zeros((NE + NI, ))
    idxvec = np.require(idxvec, np.int32, requires)
    nPostNeurons = np.require(nPostNeurons, np.int32, requires);
    # READ
    fpIdxVec = open(baseFldr + 'idxVec.dat', 'rb')
    idxvec = np.fromfile(fpIdxVec, dtype = np.int32)
    fpIdxVec.close()
    fpNpostNeurons = open(baseFldr + 'nPostNeurons.dat', 'rb')
    nPostNeurons = np.fromfile(fpNpostNeurons, dtype = np.int32)
    fpNpostNeurons.close()
    sparseVec = np.zeros((nPostNeurons.sum(), ))
    sparseVec = np.require(sparseVec, np.int32, requires)
    fpsparsevec = open(baseFldr + 'sparseConVec.dat', 'rb')
    sparseVec = np.fromfile(fpsparsevec, dtype = np.int32)
    fpsparsevec.close()
    print 'nConnections = ', nPostNeurons.sum()
    print 'sparsevec length = ', sparseVec.size


    
    nPostE = nPostNeurons[:NE].sum()
    tmp = sparseVec[sparseVec < NE]
    print np.sort(tmp[:10])

    # REWIRE
    if IF_TEST:
	po = np.linspace(0, np.pi, 1000)
    else:
	po = GetPOofPop(p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N, K, nPop, T, IF_IN_RANGE = True)
    preNeurons = Convert2InDegree(idxvec, nPostNeurons, sparseVec)
    rewiredPreNeurons = RewireSqrtK(preNeurons, p, po, K, NE, rewireType, trNo + 1)
    sparseConVec = Convert2OutDegree(rewiredPreNeurons, nPostNeurons)
    sparseVec = np.require(sparseConVec, np.int32, requires)
    tmp = sparseVec[sparseVec < NE]
    print np.sort(tmp[:10])
    print 'nConnections = ', nPostNeurons.sum()
    print 'sparsevec length = ', sparseVec.size
    # raise SystemExit
    trNo += 1 # connectivity written to another folder
    if not IF_TEST:
	baseFldrNew = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop)
    else:
	raise SystemExit
    
    os.system('mkdir -p ' + baseFldrNew)
    os.system('cp ' + baseFldr + '*.out ' + baseFldrNew)
    os.system('cp ' + baseFldr + '*FF.dat ' + baseFldrNew)

    # WRITE
    print '--'*20
    print 'saving files to: ', baseFldrNew
    fpsparsevec = open(baseFldrNew + 'sparseConVec.dat', 'wb')
    sparseVec.tofile(fpsparsevec)
    fpsparsevec.close()
    fpIdxVec = open(baseFldrNew + 'idxVec.dat', 'wb')
    idxvec.tofile(fpIdxVec)
    fpIdxVec.close()
    fpNpostNeurons = open(baseFldrNew + 'nPostNeurons.dat', 'wb')
    nPostNeurons.tofile(fpNpostNeurons)
    fpNpostNeurons.close()    
