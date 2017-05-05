import numpy as np
import pylab as plt


# def ComputeTuning(fr):
    





phis = np.arange(11.25, 180, 11.25)
phis = np.arange(22.5, 180, 22.5)
nNeurons = 1000
tc = np.zeros((nNeurons, phis.size))
for i, iPhi in enumerate(phis):
    fr = np.loadtxt('meanrates_theta%f_tr0.txt'%(iPhi))
    print fr.shape
    tc[:, i] = fr[:nNeurons]

plt.ion()
for i in np.random.randint(0, nNeurons, 101):
    plt.plot(phis, tc[i, :], 'ks-')
    plt.waitforbuttonpress()
    plt.clf()

