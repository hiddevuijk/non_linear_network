import numpy as np
from analysis import *
import sys
import warnings

warnings.filterwarnings("error")

x = np.genfromtxt("x.csv",delimiter=';')
N = x.shape[0]

xAvg = np.mean(x,axis=0)
np.savetxt('xAvg.csv',xAvg,delimiter=';')
xAvg = acorr(xAvg,norm=True)
np.savetxt('acorrAvg.csv',xAvg,delimiter=';')
meanx = np.mean(x)
x -= np.ones(x.shape)*meanx
x = acorrArray(x)
x = np.mean(x,axis=0)
np.savetxt('delta.csv',x,delimiter=';')
normalize(x)
np.savetxt('deltaN.csv',x,delimiter=';')
