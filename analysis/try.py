import numpy as np
import matplotlib.pyplot as plt
import math
from analysis import *
import sys

deltaEN = np.genfromtxt("deltaEN.csv",delimiter=';')
avgE = np.genfromtxt("averageE.csv",delimiter=';')
deltaIN = np.genfromtxt("deltaIN.csv",delimiter=';')
avgI = np.genfromtxt("averageI.csv",delimiter=';')


avgE -= np.mean(avgE)*np.ones(avgE.shape[0])
avgEAcorr= acorr(avgE,norm=True)

avgI -= np.mean(avgI)*np.ones(avgI.shape[0])
avgIAcorr = acorr(avgI,norm=True)



#plt.subplot(1,2,1)
plt.plot(avgEAcorr,color='red',label='acorr(avgE)')
plt.plot(deltaEN,color='blue',label='deltaEN')

plt.legend()
plt.xlim(0,50)
plt.ylim(-.25,1.1)
'''
plt.subplot(1,2,2)
plt.plot(avgIAcorr,color='red',label='acorr(avgI)')
plt.plot(deltaIN,color='blue',label='deltaIN')

plt.legend()
plt.xlim(0,50)
plt.ylim(-.25,1.1)
'''

plt.show()


