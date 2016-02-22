import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from shutil import copyfile

tmax = 100000

name = 'a'



what = sys.argv[1]

if len(sys.argv) >2:
	tmax = int(sys.argv[2])


ac = np.genfromtxt("../"+what+"_"+name+".csv",delimiter=';')
tt = np.genfromtxt("../t_"+name+".csv",delimiter=';')
timax = 0
while tt[timax] < tmax and timax<(len(tt) -1):
	timax += 1
plt.plot(tt[:timax],ac[:timax])
	
plt.xlabel("t")
plt.savefig(what+".pdf")

