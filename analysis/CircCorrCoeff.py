import numpy as np
from ScriptsRewire import GetPOofPop, GetPOofPopAllNeurons 

def CircularCorrCoeff(x, y):
    # x and y must be in radians
    n = x.size
    nX = x.size
    nY = y.size
    if nX != nY:
    	n = np.max([nX, nY])
    if(nX != n):
    	x = np.ones((n, )) * x
    if(nY != n):
    	y = np.ones((n, )) * y

    numerator = 0
    for i in range(n - 1):
	for j in range(i + 1, n):
	    numerator += np.sin(2.0 * (x[i] - x[j])) * np.sin(2.0 * (y[i] - y[j]))
    denom1 = 0
    denom2 = 0
    for i in range(n - 1):
	for j in range(i + 1, n):
	    denom1 += np.sin(2.0 * (x[i] - x[j]))**2
	    denom2 += np.sin(2.0 * (y[i] - y[j]))**2
    denom = np.sqrt(denom1 * denom2)
    return numerator / denom


if __name__ == "__main__":
    neuronType = 'I'
    NE = 10000
    IF_COMPUTE = False
    if IF_COMPUTE:
	if neuronType == 'E':
	    poCntr = GetPOofPop(0, 0, .075, 0.075, 'rand', kappa=0)
	    po8tr1 = GetPOofPop(0, 0, .075, 0.075, 'rand', kappa=8, trNo=1)
	    po8tr100 = GetPOofPop(0, 0, .075, 0.075, 'rand', kappa=8, trNo=100)    
	else:
	    poCntr = GetPOofPopAllNeurons(0, 0, .075, 0.075, 'rand', kappa=0)[NE:]
	    po8tr1 = GetPOofPopAllNeurons(0, 0, .075, 0.075, 'rand', kappa=8, trNo=1)[NE:]
	    po8tr100 = GetPOofPopAllNeurons(0, 0, .075, 0.075, 'rand', kappa=8, trNo=100)[NE:]   


	cccCntrVsTr1 = CircularCorrCoeff(poCntr, po8tr1)
	np.save('./data/CCC_cntrl_vs_kappa8_tr1_' + neuronType, cccCntrVsTr1)

	cccCntrVsTr100 = CircularCorrCoeff(poCntr, po8tr100)    
	np.save('./data/CCC_cntrl_vs_kappa8_tr100_' + neuronType, cccCntrVsTr100)
    else:
	


    
