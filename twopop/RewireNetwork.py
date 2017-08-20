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


# import traceback

# def info(type, value, tb):
#     traceback.print_exception(type, value, tb)
#     print
#     ipdb.pm()
# sys.excepthook = info


# from IPython import embed
# def excepthook(type, value, traceback):
#     embed()

# sys.excepthook = excepthook

# def ProcessFigure(figHdl, filepath, IF_SAVE, IF_XTICK_INT = False):
#     FixAxisLimits(figHdl)
#     axPosition = [0.25, 0.25, .65, .65]
#     paperSize = [4, 3]
#     FixAxisLimits(plt.gcf(), IF_XTICK_INT)
#     Print2Pdf(plt.gcf(), filepath, paperSize, figFormat='png', labelFontsize = 12, tickFontsize=8, titleSize = 10.0, IF_ADJUST_POSITION = True, axPosition = axPosition)
#     plt.show()


def ProcessFigure(figHdl, filepath, IF_SAVE, IF_XTICK_INT = False, figFormat = 'pdf'):
    FixAxisLimits(figHdl)
    axPosition = [0.25, 0.25, .65, .65]
    paperSize = [4, 3]
    FixAxisLimits(plt.gcf(), IF_XTICK_INT)
    Print2Pdf(plt.gcf(), filepath, paperSize, figFormat=figFormat, labelFontsize = 12, tickFontsize=8, titleSize = 10.0, IF_ADJUST_POSITION = True, axPosition = axPosition)
    plt.show()
    

def FixAxisLimits(fig, IF_XTICK_INT = False):
    ax = fig.axes[0]
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    xticks = [xmin, 0.5 * (xmin + xmax), xmax]
    if IF_XTICK_INT:
	# xticks = [int(xmin), int(0.5 *(xmin + xmax)), int(xmax)]
	ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.set_xticks(xticks)
    ax.set_yticks([ymin, 0.5 *(ymin + ymax), ymax])
    plt.draw()



#ipdb.launch_ipdb_on_exception()

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
    Ninput = len(nPostNeurons)
    Noutput = np.max(sparseConVec) + 1
    preNeurons = [[] for x in xrange(Noutput)]
    # ipdb.set_trace()
    for i in range(Ninput):
	iPostNeurons = sparseConVec[idxVec[i] : idxVec[i] + nPostNeurons[i]]
	for j, iPost in enumerate(iPostNeurons):
	    preNeurons[iPost].append(i)
    nConnections = 0
    for i in range(Noutput):
	nConnections += len(preNeurons[i])
    print '#fromPre=', nConnections, ' ', '#fromFile=', nPostNeurons.sum()
    return preNeurons

def Convert2OutDegree(preNeurons, nPostNeurons, idxvec, newLinksPost, NE):
    # sparseConVec = np.zeros((nPostNeurons.sum(), ))
    N = len(nPostNeurons) #np.max(sparseConVec) + 1    
    sparseConVec = [[] for x in xrange(N)]
    rewiredConnections = [] #[[] for x in xrange(NE)]    
    # rewiredConnections = [np.zeros((x, ), dtype=int).tolist() for x in nPostNeurons]
    count = 0
    for i in range(N):
	iPreNeurons = preNeurons[i]
	for kk in iPreNeurons:
	    sparseConVec[kk].append(i)

    for i in range(NE):
	pstNeurons = sparseConVec[i] #idxvec[i] : nPostNeurons[i] + idxvec[i]]
	tmp = np.asarray(np.in1d(pstNeurons, newLinksPost[i]), dtype = int)
	# ipdb.set_trace()
	rewiredConnections.append(tmp.tolist())
    # ipdb.set_trace()	
    return np.asarray(np.hstack(sparseConVec), dtype = np.int32), np.asarray(np.hstack(rewiredConnections), dtype = np.int32)

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
	
def RewireSqrtK(preNeurons, po, K, NE, rewireType, stepNo, kappa, nmax = 10):
    N = NE #len(preNeurons)
    recModulation = 1
    print 'nPre befor = ', len(np.hstack(preNeurons))
    linksRemoved = [[] for jjj in xrange(N)]
    newLinks = [[] for jjj in xrange(N)]
    newLinksPost = [[] for jjj in xrange(N)]   # in sparse vec representation
    for i in range(N):
	iPreNeurons = preNeurons[i]	
	iPreNeurons = np.array(iPreNeurons)
	preLenghtOld = iPreNeurons.size
	iPreNeuronsI = iPreNeurons[iPreNeurons >= N]
	iPreNeurons = iPreNeurons[iPreNeurons < N]
	iK = len(iPreNeurons)
	# ipdb.set_trace()
	if rewireType == 'rand':
            nLinks2Break = int(kappa * np.sqrt(float(iK)) / float(stepNo))	    
	    links2Break = np.sort(np.random.choice(range(iK), nLinks2Break, replace = False))
	    # print links2Break
	if rewireType == 'decay':
	    nLinks2Break = int(kappa * np.sqrt(float(iK)) / float(nmax))
	    links2Break = np.sort(np.random.choice(range(iK), nLinks2Break, replace = False))
	if rewireType == 'exp':
	    links2Break = GetNLinks(i, po, nLinks2Break, iPreNeurons)
	maxIters = nLinks2Break * 50
	if maxIters > N:
	    maxIters = N

	# vidx = np.empty((randomPreNeuron.size, ))
	# vidx[:] = True
	# ipdb.set_trace()
	randomPreNeuron = np.setdiff1d(range(N), iPreNeurons)
	randomPreNeuron = np.random.permutation(randomPreNeuron) 	
	# for mm in links2Break:
	#     vidx = np.logical_and(randomPreNeuron != iPreNeurons[mm], vidx)
	# randomPreNeuron = randomPreNeuron[vidx]
	maxIters = randomPreNeuron.size
        iterCount = 0
        nUnique = np.unique(iPreNeurons)
	linksRemoved[i] = iPreNeurons[links2Break]
	# ipdb.set_trace()
	# print '=='*50
	# print '#neuron: ', i
	# print 'preNeurons: ', iPreNeurons
	# print 'links removed: ', linksRemoved[i]
	for j in links2Break:
	    nRewiredCount = 0	
	    while nRewiredCount < 1 and iterCount < maxIters:
		# prob = 0.5 * (1 + np.cos(2.0 * (po[i] - po[randomPreNeuron[iterCount]])))
                prob = RewireProbFunc(po[i], po[randomPreNeuron[iterCount]], 'cos')		
		if(prob >= np.random.rand()):
		    # if j >= 0:
		    # if not np.any(iPreNeurons == randomPreNeuron[iterCount]):
		    newLink = randomPreNeuron[iterCount]
		    iPreNeurons[j] = newLink
		    newLinks[i].append(newLink)
		    newLinksPost[newLink].append(i)
		    nRewiredCount += 1
		iterCount += 1
        # print 'newLinks: ', newLinks[i]
	preNeurons[i] = iPreNeurons.tolist()
        # print 'rwrd preneurons ', iPreNeurons
	# print 'pst rep: ', newLinksPost[i]
	# print '=='*50
	# ipdb.set_trace()	
	preNeurons[i].extend(iPreNeuronsI.tolist())
	preLenghtNew = len(preNeurons[i])
	if(preLenghtOld != preLenghtNew):
	    print 'dose not work for neuron#', i
    print 'nPre after = ', len(np.hstack(preNeurons))    
    return preNeurons, linksRemoved, newLinks, newLinksPost

def GetBaseFolderOld(p, gamma, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, kappa = 1):
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

def GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, kappa = 1):
    # ipdb.set_trace()
    if kappa == -1:
	baseFldr = GetBaseFolderOld(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    else:
	if rewireType == 'rand':
	    tag = ''
	if rewireType == 'exp':
	    tag = '1'
	if rewireType == 'decay':
	    tag = '2'
	rootFolder = ''
	baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
	if nPop == 1:
	    baseFldr = baseFldr + 'onepop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
	#ipdb.set_trace()	
	if nPop == 2:
	    if T < 1000:
		baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), trNo)
	    else:

		if gamma >= .1 or gamma == 0:
		    baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
		else:
		    baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(T*1e-3), trNo)
    return baseFldr

# def GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, kappa = 1):
#     # ipdb.set_trace()
#     if rewireType == 'rand':
# 	tag = ''
#     if rewireType == 'exp':
# 	tag = '1'
#     rootFolder = ''
#     baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
#     if nPop == 1:
#     	baseFldr = baseFldr + 'onepop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
#     # ipdb.set_trace()	
#     if nPop == 2:
#         if T < 1000:
# 	    baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), trNo)
#         else:
# 	    if gamma >= .1 or gamma == 0:
# 		baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
# 	    else:
# 		baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(T*1e-3), trNo)
#     return baseFldr

# def GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, kappa = 1):
#     if rewireType == 'rand':
# 	tag = ''
#     if rewireType == 'exp':
# 	tag = '1'
#     rootFolder = ''
#     baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
#     if nPop == 1:
#     	baseFldr = baseFldr + 'onepop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
#     if nPop == 2:
# 	if gamma >= .1 or gamma == 0:
# 	    baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
# 	else:
# 	    baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/p%sgamma/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(T*1e-3), trNo)
#     print baseFldr	    
#     return baseFldr
		
def LoadFr(p, gamma, phi, mExt, mExtOne, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False, rewireType = 'rand', kappa = 1):
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    if IF_VERBOSE:
    	print baseFldr
    filename = 'meanrates_theta%.6f_tr%s.txt'%(phi, trNo)
    print filename
    return np.loadtxt(baseFldr + filename)

def GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 1):
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
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, trNo, T, NE, K, nPop, rewireType = rewireType, IF_VERBOSE = True, kappa = kappa)
	    else:
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, trNo, T, NE, K, nPop, rewireType = rewireType, IF_VERBOSE = False, kappa = kappa)
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

def GetPOofPop(p, gamma, mExt, mExtOne, rewireType, nPhis = 8, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_IN_RANGE = True, kappa = 1):
    nNeurons = N
    thetas= np.linspace(0, np.pi, nNeurons, endpoint = False)
    tc = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
    prefferedOri = POofPopulation(tc[:N], IF_IN_RANGE = True) * np.pi / 180.0
    return prefferedOri

def ComputeFFInput_OLD(nPhis, kappa, p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, NFF = 10000, JE0 = 2.0, JI0 = 1.0):
    idxvecFF, nPostNeuronsFF, sparseVecFF = LoadSparseMatFF(kappa, p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop)
    NE = N
    NI = N
    inputFF = np.zeros((NE + NI, ))
    preNeuronsFF = Convert2InDegree(idxvecFF, nPostNeuronsFF, sparseVecFF)
    uofPhi = np.zeros((NE + NI, nPhis))
    thetas = np.linspace(0, np.pi, nPhis, endpoint = False)
    ffIdx = np.arange(NFF) * np.pi / NFF
    mFF = np.zeros((NFF, ))
    # ipdb.set_trace()
    for k, thetaExt in enumerate(thetas):
	for ffi in range(NFF):
	    mFF[ffi] = mExt + mExtOne * np.cos(2.0 * (thetaExt - ffIdx[ffi]))
        ipdb.set_trace()
	# for i in range(NE):
	#     for j in 
	#     inputFF[i] = 
	    
	# for i in xrange(NE):
	#     inputFF[i] = JE0 * np.sum(mFF[preNeuronsFF[i]])
	for i in range(NE, NE + NI):
	    inputFF[i] = JI0 * np.sum(mFF[preNeuronsFF[i]]) 
	uofPhi[:, k] = inputFF
    return uofPhi

def ComputeFFInput(nPhis, p = 0, gamma = 0, kappa = 1, mExt = .075, mExtOne= .075, trNo=0, N=10000, K=1000, nPop=2, NFF = 10000, JE0 = 2.0, JI0 = 1.0, cFF = 0.2, rewireType = 'rand', T = 1000):
    print 'in rw'
    idxvecFF, nPostNeuronsFF, sparseVecFF = LoadSparseMatFF(kappa, p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop) 
    NE = N
    NI = N
    JE0 = JE0 / (cFF * np.sqrt(K))
    JI0 = JI0 / (cFF * np.sqrt(K))
    inputFF = np.zeros((NE + NI, ))
    preNeuronsFF = Convert2InDegree(idxvecFF, nPostNeuronsFF, sparseVecFF)
    uofPhi = np.zeros((NE + NI, nPhis))
    thetas = np.linspace(0, np.pi, nPhis, endpoint = False)
    ffIdx = np.arange(NFF) * np.pi / NFF
    mFF = np.zeros((NFF, ))
    # ipdb.set_trace()
    for k, thetaExt in enumerate(thetas):
	for ffi in range(NFF):
	    mFF[ffi] = mExt + mExtOne * np.cos(2.0 * (thetaExt - ffIdx[ffi]))
        # ipdb.set_trace()
	for i in xrange(NE):
	    inputFF[i] = JE0 * np.sum(mFF[preNeuronsFF[i]])
	for i in range(NE, NE + NI):
	    inputFF[i] = JI0 * np.sum(mFF[preNeuronsFF[i]]) 
	uofPhi[:, k] = inputFF
    return uofPhi

def LoadSparseMatFF(kappa, p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop):
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    print 'in rw', baseFldr
    requires = ['CONTIGUOUS', 'ALIGNED']
    # READ
    fpIdxVec = open(baseFldr + 'idxVecFF.dat', 'rb')
    idxvec = np.fromfile(fpIdxVec, dtype = np.int32)
    fpIdxVec.close()
    fpNpostNeurons = open(baseFldr + 'nPostNeuronsFF.dat', 'rb')
    nPostNeurons = np.fromfile(fpNpostNeurons, dtype = np.int32)
    fpNpostNeurons.close()
    sparseVec = np.zeros((nPostNeurons.sum(), ))
    sparseVec = np.require(sparseVec, np.int32, requires)
    fpsparsevec = open(baseFldr + 'sparseConVecFF.dat', 'rb')
    sparseVec = np.fromfile(fpsparsevec, dtype = np.int32)
    fpsparsevec.close()
    return idxvec, nPostNeurons, sparseVec
        
def LoadSparseMat(kappa, p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop):
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
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

def ConTest(n = 1000, K = 100):
    x = np.pi * np.random.rand(n)
    y = np.pi * np.random.rand(n)
    z = np.zeros((n, ))
    poDiffTmp = []
    Ki = np.floor(np.random.randn(n) * np.sqrt(K) + K)
    # ipdb.set_trace()
    for i in range(n):
	preNeurons = np.random.choice(n, Ki[i], replace=False)
	z[i] = np.abs(np.sum(np.exp(2j * (x[i] - y[preNeurons])))) / float(Ki[i])
	poDiffTmp.append(AngleDiff(x[i], y[preNeurons]))
    # plt.hist(z, 100, normed = True, histtype = 'step', label = 'N = %s, K = %s'%(n, K))
    poDiff = []
    for i in range(n):
	for k in poDiffTmp[i]:
	    poDiff.append(k)
    plt.hist(poDiff, 100, normed = 1, histtype = 'step')
    return poDiff
    
def ConnectionDistr(kappa, p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N = 10000, K = 1000, nPop = 2, T = 1000):
    po = GetPOofPop(p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N, K, nPop, T, IF_IN_RANGE = True, kappa = kappa)
    idxvec, nPostNeurons, sparseVec = LoadSparseMat(kappa, p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop)
    preNeurons = Convert2InDegree(idxvec, nPostNeurons, sparseVec)
    z = np.empty((N, ))
    z[:] = np.nan
    poDiff = []
    poDiffCirc = []
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
		poDiff.append(poi - pok)
		poDiffCirc.append(AngleDiff(poi, pok))
	z[i] = np.abs(tmp) / float(Ki)
    return z, poDiff, poDiffCirc

def AngleDiff(v1, v2, IF_IN_DEGREES = False):
    if IF_IN_DEGREES:
	v1 = np.pi * v1 / 180.0
	v2 = np.pi * v2 / 180.0
    v1x = np.cos(2 * v1);
    v1y = np.sin(2 * v1);
    v2x = np.cos(2 * v2);
    v2y = np.sin(2 * v2);
    out = np.arccos(v1x * v2x + v1y* v2y) * 0.5;
    if IF_IN_DEGREES:
	out = out * 180 / np.pi
    return out

def CompareConDistr(trList, kappa = 1, p = 0, gamma = 0, mExt = 0.075, mExtOne = 0.075, rewireType = 'rand', nPhis = 8, N = 10000, K = 1000, nPop = 2, T = 1000):
 #   fg0, ax0 = plt.subplots()
  #  fg1, ax1 = plt.subplots()
    fg2, ax2 = plt.subplots()
    plt.ion()    
    for trNo in trList:
	z, poDiff, poDiffCirc = ConnectionDistr(kappa, p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N, K, nPop, T)
#	ax0.hist(z, 50, normed = True, histtype = 'step', label = 'STP: %s'%(trNo))
	# ipdb.set_trace()
#	ax1.hist(np.asarray(poDiff) * 180 / np.pi, 100, normed = True, histtype = 'step', label = 'STP: %s'%(trNo))
	poDiffCirc = np.asarray(poDiffCirc)
	ax2.hist(poDiffCirc[~np.isnan(poDiffCirc)] * 180 / np.pi, 100, normed = True, histtype = 'step', label = 'STP: %s'%(trNo))
	plt.show()
    # ax0.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 8})
    # ax0.set_ylabel('Density')
    # ax0.set_xlabel(r"$\langle \exp (\phi^{po}_i - \phi^{po}_j) \rangle_j$")
    # ax1.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 8})
    # ax1.set_ylabel('Density')
    # ax1.set_xlabel(r"$\Delta PO$")

    ax2.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 8})
    ax2.set_ylabel('Density')
    ax2.set_xlabel(r"$\Delta PO \mathrm{(deg)}$")
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    # filepath = './figs/conDistr_p%sg%s'%(int(p * 100), int(gamma * 100))
    # ProcessFigure(fg0, filepath, True)
    # filepath = './figs/conDistr_vs_deltaPOs_p%sg%s'%(int(p * 100), int(gamma * 100))
    # ProcessFigure(fg1, filepath, True)
    filepath = './figs/conDistr_DeltaPOCirc_p%sg%s_'%(int(p * 100), int(gamma * 100)) + rewireType
    ProcessFigure(fg2, filepath, True)
    ipdb.set_trace()



def GetPODiff(kappa, p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N = 10000, K = 1000, nPop = 2, T = 1000):
    po = GetPOofPop(p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N, K, nPop, T, IF_IN_RANGE = True, kappa = kappa)
    idxvec, nPostNeurons, sparseVec = LoadSparseMat(kappa, p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop)
    preNeurons = Convert2InDegree(idxvec, nPostNeurons, sparseVec)
    out = []
    for i in range(N):
	preNeuronsTmp = np.asarray(preNeurons[i], dtype = np.int)
	preNeuronsTmp = preNeuronsTmp[preNeuronsTmp < N]
	for k in preNeuronsTmp:
	    # out.append(np.abs(po[i] - po[k]))
	    out.append(np.abs(po[k]))	    
    return np.asarray(out, dtype = np.float)

def CompareDPOvsDist(trList, kappa = 1, p = 0, gamma = 0, mExt = 0.075, mExtOne = 0.075, rewireType = 'rand', nPhis = 8, N = 10000, K = 1000, nPop = 2, T = 1000):
    plt.figure()
    plt.ion()
    for trNo in trList:
	z = GetPODiff(kappa, p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N, K, nPop, T)
	plt.hist(z, 50, normed = True, histtype = 'step', label = 'STP: %s'%(trNo))
	plt.show()
    plt.legend(loc = 0, frameon = False, numpoints = 1, prop = {'size': 8})
    plt.ylabel('Density')
    plt.xlabel(r"$\langle \phi^{po}_i - \phi^{po}_j \rangle$")
    plt.savefig('./figs/DeltaPO_vs_conDistr_p%sg%s.png'%(int(p * 100), int(gamma * 100)))
  
if __name__ == '__main__':
    # NetworkType : {'uni', 'ori'}, 'uni' is for standard random network, 'ori' is to rewire depending on the distance in ori space
    [trNo, kappa, K, rewireType, p, gamma, mExt, mExtOne, nPhis,NE, NI, nPop, T] = DefaultArgs(sys.argv[1:], [0, 1, 500, 'rand', 0, 0, .075, .075, 8, 10000, 10000,  2, 1000])
    NE = int(NE)
    NI = int(NI)
    N = NE
    K = int(K)
    p = float(p)
    kappa = float(kappa)
    trNo = int(trNo)
    mExt = float(mExt)
    mExtOne = float(mExtOne)
    gamma = float(gamma)
    nPop = int(nPop)
    nPhis = int(nPhis)
    T = int(T)
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    IF_TEST = False
    if IF_TEST:
	baseFldr = './'
	NE = 10 
	NI = 10
	K = 2
    requires = ['CONTIGUOUS', 'ALIGNED']
    # convec = np.require(convec, np.int32, requires)
    nPostNeurons = np.zeros((NE + NI, ))
    idxvec = np.zeros((NE + NI, ))
    idxvec = np.require(idxvec, np.int32, requires)
    nPostNeurons = np.require(nPostNeurons, np.int32, requires)
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
	po = np.linspace(0, np.pi, NE)
    else:
	po = GetPOofPop(p, gamma, mExt, mExtOne, rewireType, nPhis, trNo, N, K, nPop, T, IF_IN_RANGE = True, kappa = kappa)
    preNeurons = Convert2InDegree(idxvec, nPostNeurons, sparseVec)
    # print preNeurons[9]
    # ipdb.set_trace()
    rewiredPreNeurons, linksRemoved, newLinks, newLinksPost  = RewireSqrtK(preNeurons, po, K, NE, rewireType, trNo + 1, kappa)
    sparseConVec, newPostNeurons = Convert2OutDegree(rewiredPreNeurons, nPostNeurons, idxvec, newLinksPost, NE)
    # ipdb.set_trace()        
    sparseVec = np.require(sparseConVec, np.int32, requires)
    tmp = sparseVec[sparseVec < NE]
    print np.sort(tmp[:10])
    print 'nConnections = ', nPostNeurons.sum()
    print 'sparsevec length = ', sparseVec.size
    # raise SystemExit

    np.save(baseFldr + 'newLinks', newLinks)
    np.save(baseFldr + 'linksRemoved', linksRemoved)
    
    trNo += 1 # connectivity written to another folder
    if not IF_TEST:
	baseFldrNew = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    else:
	raise SystemExit
    
    os.system('mkdir -p ' + baseFldrNew)
    os.system('cp ' + baseFldr + '*.out ' + baseFldrNew)
    os.system('cp ' + baseFldr + '*FF.dat ' + baseFldrNew)

    # WRITE
    print '--'*20
    print 'saving files to: ', baseFldrNew
    print '--'*20
    print '# new connctions', newPostNeurons.sum(), '# EI_IE cons=', nPostNeurons[:NE].sum()
    fpNewPost = open(baseFldrNew + 'newPostNeurons.dat', 'wb')
    newPostNeurons = np.require(newPostNeurons, np.int32, requires)    
    newPostNeurons.tofile(fpNewPost)
    fpNewPost.close()
    fpsparsevec = open(baseFldrNew + 'sparseConVec.dat', 'wb')
    sparseVec.tofile(fpsparsevec)
    fpsparsevec.close()
    fpIdxVec = open(baseFldrNew + 'idxVec.dat', 'wb')
    idxvec.tofile(fpIdxVec)
    fpIdxVec.close()
    fpNpostNeurons = open(baseFldrNew + 'nPostNeurons.dat', 'wb')
    nPostNeurons.tofile(fpNpostNeurons)
    fpNpostNeurons.close()    

