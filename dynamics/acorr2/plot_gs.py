import matplotlib.pyplot as plt
import numpy as np
import sys


tplot = 50
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

#	tick = np.linspace(0,1,tplot)*t[tplot -1]
#	plt.plot(tick,ac[:tplot],label=g)
	tt = 0
	while t[tt] < tplot:
		tt+=1
	tick = np.linspace(0,1,tt)*t[tt]
	plt.plot(tick,ac[:tt],label=g)

plt.legend()
if(not save):
	plt.show()
else:
	plt.savefig(name)




