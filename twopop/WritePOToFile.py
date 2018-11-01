import numpy as np
import sys

def GetBaseFolder_OLD(p, gamma, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, kappa = 1):
    # ipdb.set_trace()
    tag = ''
    if kappa == -1:
	baseFldr = GetBaseFolderOld(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    else:
	if rewireType == 'rand':
	    tag = ''
	if rewireType == 'exp':
	    tag = '1'
	if rewireType == 'decay':
	    tag = '2'
	if rewireType == 'twosteps':
	    tag = '3'	    
	rootFolder = ''
	baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
	if nPop == 1:
	    baseFldr = baseFldr + 'onepop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
	# ipdb.set_trace()	
	if nPop == 2:
	    if T < 1000:
		baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), trNo)
	    else:
		if gamma >= .1 or gamma == 0:
		    baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
		else:
		    baseFldr = baseFldr + 'twopop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(T*1e-3), trNo)
    if rewireType == 'cntrl':
        baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
        baseFldr = baseFldr + 'twopop/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    return baseFldr


def GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, kappa = 1):
    # ipdb.set_trace()
    tag = ''
    if kappa == -1:
	baseFldr = GetBaseFolderOld(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    else:
	if rewireType == 'rand':
	    tag = ''
	if rewireType == 'exp':
	    tag = '1'
	if rewireType == 'decay':
	    tag = '2'
	if rewireType == 'twosteps':
	    tag = '3'	    
	rootFolder = ''
	baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
	if nPop == 1:
	    baseFldr = baseFldr + 'onepop/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
	# ipdb.set_trace()	
	if nPop == 2:
	    if T < 1000:
		baseFldr = baseFldr + 'twopop/PRX_data/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), trNo)
	    else:
		if gamma >= .1 or gamma == 0:
		    baseFldr = baseFldr + 'twopop/PRX_data/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma%s/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
		else:
		    baseFldr = baseFldr + 'twopop/PRX_data/data/rewire%s/N%sK%s/m0%s/mExtOne%s/kappa%s/p%sgamma/T%s/tr%s/'%(tag, N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(kappa * 10), int(p * 10), int(T*1e-3), trNo)
    if rewireType == 'cntrl':
        baseFldr = rootFolder + '/homecentral/srao/Documents/code/binary/c/'
        baseFldr = baseFldr + 'twopop/PRX_data/data/N%sK%s/m0%s/mExtOne%s/p%sgamma%s/T%s/tr%s/'%(N, K, int(1e3 * mExt), int(1e3 * mExtOne), int(p * 10), int(gamma * 10), int(T*1e-3), trNo)
    return baseFldr
    

def LoadFr(p, gamma, phi, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False, kappa = 0):
    #ipdb.set_trace()
    baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    filename = 'meanrates_theta%.6f_tr%s.txt'%(phi, trNo)
    if IF_VERBOSE:
    	print baseFldr
	print filename
    return np.loadtxt(baseFldr + filename)


def GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 0, IF_SUCCESS = False):
    NE = N
    NI = N
    tc = np.zeros((NE + NI, nPhis))
    tc[:] = np.nan
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    IF_FILE_LOADED = False
    for i, iPhi in enumerate(phis):
	print i, iPhi
	try:
	    if i == 0:
		print 'loading from fldr: ',	    
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = True, kappa = kappa)
	    else:
		fr = LoadFr(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = False, kappa = kappa)
	    if(len(fr) == 1):
		if(np.isnan(fr)):
		    print 'file not found!'
	    tc[:, i] = fr
            IF_FILE_LOADED = True
	except IOError:
	    print 'file not found!'
	    # raise SystemExit
    if IF_SUCCESS:
        return tc, IF_FILE_LOADED
    else:
        return tc


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
    zk = np.dot(firingRate, np.exp(2j * atTheta * np.pi / 180))
    out = np.angle(zk) * 180.0 / np.pi
    if IF_IN_RANGE:
	if(out < 0):
	    out += 360
    return out * 0.5

        
def WritePOToTr0FolderBeforeRewiring(kappa, mExt, mExtOne, trNo, p = 0, gamma = 0, nPhis=8, rewireType='rand', N=10000, K=1000, nPop=2, T=1000):
    out = np.nan
    try:
	tcOut = GetTuningCurves(p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
        po = POofPopulation(tcOut, IF_IN_RANGE = 1)
        requires = ['CONTIGUOUS', 'ALIGNED']
        po = np.require(po, np.float64, requires) * np.pi / 180.0
        baseFldr = GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
        fp = open(baseFldr + 'poOfNeurons.dat', 'wb')
        po.tofile(fp)
        fp.close()
        print po[:11]
        print 'file writen to: ', baseFldr
    except IOError:
        print 'file not found'

if __name__ == "__main__":
    kappa = float(sys.argv[1])
    mExt = float(sys.argv[2])
    mExtOne = float(sys.argv[3])
    trNo = int(sys.argv[4])
    K=int(sys.argv[5])
    T=int(sys.argv[6])
    WritePOToTr0FolderBeforeRewiring(kappa, mExt, mExtOne, trNo, K=K, T=T)
    
    
    
