import matplotlib.pyplot as plt
import numpy as np
import sys

save = False
name=''
if len(sys.argv) == 2:
	name = sys.argv[1]
	save = True

names = np.loadtxt("names.txt",dtype=np.str)
for n in names:
#	xt = np.genfromtxt("x_" + n + ".csv",delimiter=';')
#	t = np.genfromtxt("t_" + n + ".csv",delimiter=';')
	ac = np.genfromtxt("acorr_"+n+".csv",delimiter=';')
	acm = np.genfromtxt("acorr_mean_"+n+".csv",delimiter=';')

	plt.subplot(1,2,1)
	plt.plot(ac,label=n)
	plt.title("mean of autocorr")

	plt.subplot(1,2,2)
	plt.plot(acm,label=n)
	plt.title("autocorr of mean")

plt.legend()
if(not save):
	plt.show()
else:
	plt.savefig(name)




