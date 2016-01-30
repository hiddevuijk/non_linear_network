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
if len(sys.argv) == 2:
	name = sys.argv[1]
	save = True

gs = np.loadtxt("names.txt",dtype=np.str)

for n in range(len(gs)):
	xt = np.genfromtxt("x_" + gs[n] +".csv",delimiter=';')
	t = np.genfromtxt("t_"+gs[n]+".csv",delimiter=';')

	tf = t[-1]
	N= xt.shape[0]
	tmax=xt.shape[1]
	dt = tf/tmax

	ac = np.asarray([0.0]*((tmax-100)/2))
	for i in range(N):
		ac += acorr(xt[i][100:])/N
	N = min(N,3)

	plt.subplot(len(gs),1,n+1)
	plt.title(gs[n])
	for i in range(N):
		plt.plot(t,xt[i])

if(not save):
	plt.show()
else:
	plt.savefig(name)




