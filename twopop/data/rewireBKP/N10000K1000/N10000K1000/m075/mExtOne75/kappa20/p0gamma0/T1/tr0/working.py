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
