import matplotlib.pyplot as plt
import numpy as np
import sys


def acorr(a):
	a -= np.mean(a)
	ac = np.correlate(a,a,"same")[len(a)/2:]
	ac0 = ac[0]
	return ac/ac0



save = False
name=''
name_out=''
if len(sys.argv) == 2:
	name = "_" + sys.argv[1]

if len(sys.argv) == 3:
	name_out = sys.argv[2]
	save = True

xt = np.genfromtxt("x"+name+".csv",delimiter=';')
t = np.genfromtxt("t"+name+".csv",delimiter=';')

tf = t[-1]
N= xt.shape[0]


N = min(N,3)

for i in range(N):
	plt.plot(t,xt[i])



if(not save):
	plt.show()
else:
	plt.svefig(name_out)




