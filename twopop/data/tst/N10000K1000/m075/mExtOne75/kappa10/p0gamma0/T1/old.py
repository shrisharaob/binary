def ComputeFFInput(nPhis, kappa, p, gamma, mExt, mExtOne, rewireType, trNo, T, N, K, nPop, NFF = 10000, JE0 = 2.0, JI0 = 1.0):
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
        # ipdb.set_trace()
	for i in xrange(NE):
	    inputFF[i] = JE0 * np.sum(mFF[preNeuronsFF[i]])
	for i in range(NE, NE + NI):
	    inputFF[i] = JI0 * np.sum(mFF[preNeuronsFF[i]]) 
	uofPhi[:, k] = inputFF
    return uofPhi
