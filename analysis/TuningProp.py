import numpy as np
import pylab as plt


def LoadFr(p, gamma, phi, N = 10000, K = 1000, nPop = 2):
    baseFldr = '/homecentral/srao/Documents/code/binary/c/'
    if nPop == 1:
	baseFldr = baseFldr + 'onepop/data/N%sK%s/'%(N, K)
    if nPop == 2:
	baseFldr = baseFldr + 'twopop/data/N%sK%s/'%(N, K)

    return np.loadtxt(baseFldr + 'meanrates_theta%f_tr0.txt'%(phi))

def ComputeTuning(p, gamma, nPhis, N = 10000, K = 1000, nPop = 2):
    NE = N
    NI = N
    tc = np.zeros((NE + NI, nPhis))
    phis = np.linspace(0, 180, nPhis, endpoint = False)
    for i, iPhi in enumerate(phis):
	fr = LoadFr(p, gamma, iPhi, NE, K, nPop)
	print fr.shape
	tc[:, i] = fr[:nNeurons]
    plt.ion()
    for i in np.random.randint(0, 10000, 101):
	plt.plot(phis, tc[i, :], 'ks-')
	plt.waitforbuttonpress()
	plt.clf()



    
