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

def LoadFr(phi, trNo = 0,  IF_VERBOSE = 1):# p, gamma, phi, mExt, mExtOne, rewireType, trNo = 0, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False, kappa = 1):
    # ipdb.set_trace()
    baseFldr = '' #GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
    if IF_VERBOSE:
    	print baseFldr
    filename = 'meanrates_theta%.6f_tr%s.txt'%(phi, trNo)
    print filename
    return np.loadtxt(baseFldr + filename)



def GetTuningCurves(p = 0, gamma = 0, nPhis = 8, mExt = .075, mExtOne=.075, rewireType='rand', trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 1):
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
		fr = LoadFr(iPhi, trNo) #p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = True, kappa = kappa)
	    else:
		fr = LoadFr(iPhi, trNo) #(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = False, kappa = kappa)
	    if(len(fr) == 1):
		if(np.isnan(fr)):
		    print 'file not found!'
	    tc[:, i] = fr
	except IOError:
	    print 'file not found!'
	    raise SystemExit
    return tc

def POofPopulation(tc, theta = np.arange(0.0, 180.0, 22.5), IF_IN_RANGE = False):
    # return value in degrees
    nNeurons, nAngles = tc.shape
    theta = np.linspace(0, 180, nAngles, endpoint = False)
    po = np.zeros((nNeurons, ))
    for kNeuron in np.arange(nNeurons):
        po[kNeuron] = GetPhase(tc[kNeuron, :], theta, IF_IN_RANGE)
    return po

def GetPOofPop(nPhis = 8, trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, IF_IN_RANGE = True, kappa = 1):
    nNeurons = N
    thetas= np.linspace(0, np.pi, nNeurons, endpoint = False)
    tc = GetTuningCurves() #p, gamma, nPhis, mExt, mExtOne, rewireType, trNo, N, K, nPop, T, kappa)
    prefferedOri = POofPopulation(tc[:N], IF_IN_RANGE = True) * np.pi / 180.0
    return prefferedOri

def GetPhase(firingRate, atTheta, IF_IN_RANGE = False):
    out = np.nan
    zk = np.dot(firingRate, np.exp(2j * atTheta * np.pi / 180))
    out = np.angle(zk) * 180.0 / np.pi
    if IF_IN_RANGE:
	if(out < 0):
	    out += 360
    return out * 0.5

def GetInputTuningCurves(p=0, gamma = 0, nPhis = 8, mExt = .075, mExtOne = .075, rewireType = 'rand', trNo = 0, N = 10000, K = 1000, nPop = 2, T = 1000, kappa = 1):
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
		fr = LoadInput(iPhi, trNo) #p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = True, kappa = kappa)
	    else:
		fr = LoadInput(iPhi, trNo) #LoadInput(p, gamma, iPhi, mExt, mExtOne, rewireType, trNo, T, NE, K, nPop, IF_VERBOSE = False, kappa = kappa)
	    if(len(fr) == 1):
		if(np.isnan(fr)):
		    print 'file not found!'
	    tc[:, i] = fr
	except IOError:
	    print 'file not found!'
	    raise SystemExit
    return tc


def LoadInput(phi, trNo = 0,  IF_VERBOSE = 1): #, T = 1000, N = 10000, K = 1000, nPop = 2, IF_VERBOSE = False, kappa = 1):
    # ipdb.set_trace()
    baseFldr = '' #GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop)
    if IF_VERBOSE:
    	print baseFldr
    filename = 'meanFFinput_theta%.6f_tr%s_last.txt'%(phi, trNo)
    print filename
    # ipdb.set_trace()
    return np.loadtxt(baseFldr + filename)


def ComputeFFInputNew(nPhis, mExt = .075, mExtOne= .075, trNo=0, N=10000, K=1000, nPop=2, NFF = 10000, JE0 = 2.0, JI0 = 1.0):
    idxvecFF, nPostNeuronsFF, sparseVecFF = LoadSparseMatFF()
    NE = N
    NI = N
    JE0 = JE0 / (.2 * K)
    JI0 = JI0 / (.2 * K)    
    inputFF = np.zeros((NE + NI, ))
    uofPhi = np.zeros((NE + NI, nPhis))
    thetas = np.linspace(0, np.pi, nPhis, endpoint = False)
    ffIdx = np.arange(NFF) * np.pi / NFF
    mFF = np.zeros((NFF, ))
    # ipdb.set_trace()
    for k, thetaExt in enumerate(thetas):
	for ffi in range(NFF):
	    mFF[ffi] = mExt + mExtOne * np.cos(2.0 * (thetaExt - ffIdx[ffi]))
        # ipdb.set_trace()
	for i in xrange(NFF):
	    postNeurons = sparseVecFF[idxvecFF[i] : idxvecFF[i] + nPostNeuronsFF[i]]
	    postE = postNeurons[postNeurons < N]
	    for l in postE:
		inputFF[l] += JE0
	uofPhi[:, k] = inputFF
        inputFF = np.zeros((NE + NI, ))
    return uofPhi

def ComputeFFInput(nPhis, mExt = .075, mExtOne= .075, trNo=0, N=10000, K=1000, nPop=2, NFF = 10000, JE0 = 2.0, JI0 = 1.0):
    idxvecFF, nPostNeuronsFF, sparseVecFF = LoadSparseMatFF() 
    NE = N
    NI = N
    JE0 = JE0 / (0.2 * np.sqrt(K))
    JI0 = JI0 / (0.2 * np.sqrt(K))
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


def LoadSparseMatFF(trNo = 0, N = 10000):
    baseFldr = '' #GetBaseFolder(p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, kappa)
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
