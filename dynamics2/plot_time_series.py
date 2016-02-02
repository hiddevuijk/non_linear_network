import matplotlib.pyplot as plt
import numpy as np
import sys


def acorr(a):
	a -= np.mean(a)
	ac = np.correlate(a,a,"same")[len(a)/2:]
	ac0 = ac[0]
	return ac/ac0

def psd(a):
	a = np.fft.fft(a)
	a=a[:len(a)/2]
	return abs(a)



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

plt.subplot(1,3,1)
for i in range(N):
	plt.plot(t,xt[i])

plt.subplot(1,3,2)
plt.plot(abs(np.fft.fft(acorr(xt[1]))))

plt.subplot(1,3,3)
plt.plot(psd(xt[1]))


if(not save):
	plt.show()
else:
	plt.svefig(name_out)




