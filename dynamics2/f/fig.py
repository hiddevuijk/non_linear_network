import matplotlib.pyplot as plt
import numpy as np
import sys

save = False
name=''
if len(sys.argv) == 2:
	name = sys.argv[1]
	save = True

def acorr(a):
	a -= np.mean(a)
	ac = np.correlate(a,a,"same")[len(a)/2:]
	ac0 = ac[0]
	return ac/ac0

def f(mE):
	return 1./(mE+1)

names = np.loadtxt("names.txt",dtype=np.str)

'''
n = names[0]
ac = np.genfromtxt("acorr_"+n+".csv",delimiter=';')
xt = np.genfromtxt("x_" + n + ".csv",delimiter=';')
N = xt.shape[0] 

c = np.asarray([0.0]*(xt.shape[1]/2))
for i in range(N):
	c += acorr(xt[i])/N
'''

for n in names:
#	xt = np.genfromtxt("x_" + n + ".csv",delimiter=';')
#	t = np.genfromtxt("t_" + n + ".csv",delimiter=';')
	ac = np.genfromtxt("acorr_"+n+".csv",delimiter=';')
	acm = np.genfromtxt("acorr_mean_"+n+".csv",delimiter=';')

	
	meanE = float(n)
	ff = f(meanE)
	plt.subplot(1,2,1)
	plt.plot(ac,label=ff)
	plt.title("mean of autocorr")

	plt.subplot(1,2,2)
	plt.plot(acm,label=ff)
	plt.title("autocorr of mean")

plt.legend()
if(not save):
	plt.show()
else:
	plt.savefig(name)




