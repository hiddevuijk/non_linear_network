import matplotlib.pyplot as plt
import numpy as np
import sys

save = False
name=''
if len(sys.argv) == 2:
	name = sys.argv[1]
	save = True



names = np.loadtxt("names.txt",dtype=np.str)
n = names[0]

for n in names:
	psd = np.genfromtxt("psd_"+n+".csv")
	f = np.genfromtxt("freq_"+n+".csv")
	plt.plot(psd,label=n)

plt.legend()
if(not save):
	plt.show()
else:
	plt.savefig(name)




