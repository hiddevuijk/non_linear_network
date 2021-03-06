import matplotlib.pyplot as plt
import numpy as np
import sys



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
ac = np.genfromtxt("acorrn"+name+".csv",delimiter=';')
acm = np.genfromtxt("acorr_mean"+name+".csv",delimiter=';')
m=np.genfromtxt("mean"+name+".csv",delimiter=';')

plt.plot(ac)



if(not save):
	plt.show()
else:
	plt.svefig(name_out)




