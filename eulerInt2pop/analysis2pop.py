import numpy as np
from analysis import *
import sys
import warnings

warnings.filterwarnings("error")

Ne = 0
Ni = 0

tiMax = 0
try:
	xE = np.genfromtxt("xE.csv",delimiter=';')
	Ne = xE.shape[0]
	tiMax = xE.shape[1]
except:
	pass
try:
	xI = np.genfromtxt("xI.csv",delimiter=';')
	Ni = xI.shape[0]
	tiMax = xI.shape[1]
except:
	pass


if Ne >0:
	xEAvg = np.mean(xE,axis=0)
	np.savetxt('xEAvg.csv',xEAvg,delimiter=';')
	xEAvg = acorr(xEAvg,norm=True)
	np.savetxt('acorrEAvg.csv',xEAvg,delimiter=';')

	meanxE = np.mean(xE)
	xE -= np.ones(xE.shape)*meanxE
	#xE -= np.asarray([np.mean(xE,axis=1)]*xE.shape[1]).T
	xE = acorrArray(xE)
	xE = np.mean(xE,axis=0)
	np.savetxt('deltaE.csv',xE,delimiter=';')
	normalize(xE)
	np.savetxt('deltaEN.csv',xE,delimiter=';')

if Ni>0:
	xIAvg = np.mean(xI,axis=0)
	np.savetxt('xIAvg.csv',xIAvg,delimiter=';')
	xIAvg = acorr(xIAvg,norm=True)
	np.savetxt('acorrIAvg.csv',xIAvg,delimiter=';')
	#meanxI = np.mean(xI)
	#xI -= np.ones(xI.shape)*meanxI
	xI -= np.asarray([np.mean(xI,axis=1)]*xI.shape[1]).T
	xI = acorrArray(xI)
	xI = np.mean(xI,axis=0)
	np.savetxt('deltaI.csv',xI,delimiter=';')
	normalize(xI)
	np.savetxt('deltaIN.csv',xI,delimiter=';')



