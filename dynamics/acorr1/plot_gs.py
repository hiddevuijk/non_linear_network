import matplotlib.pyplot as plt
import numpy as np
import sys

def acorr(a):
	a -= np.mean(a)
	ac = np.correlate(a,a,"same")[len(a)/2:]
	return ac/ac[0]

save = False
name=''
if len(sys.argv) == 2:
	name = sys.argv[1]
	save = True

gs = np.loadtxt("names.txt",dtype=np.str)
for g in gs:
	xt = np.genfromtxt("x_" + g + ".csv",delimiter=';')
	t = np.genfromtxt("t_" + g + ".csv",delimiter=';')

	tf = t[-1]
	N= xt.shape[0]
	tmax=xt.shape[1]
	dt = tf/tmax

	ac = np.asarray([0.0]*((tmax-100)/2))
	for i in range(N):
		ac =ac + acorr(xt[i][100:])/N

	tick = np.linspace(0,1,tmax/4)*t[tmax/4 -1]
	plt.plot(tick,ac[:tmax/4],label=g)
plt.legend()
if(not save):
	plt.show()
else:
	plt.savefig(name)




