import matplotlib.pyplot as plt
import numpy as np
import sys

save = False
name=''
if len(sys.argv) == 2:
	name = sys.argv[1]
	save = True


r = 15

re = np.genfromtxt("eval_re.csv",delimiter=';')
im = np.genfromtxt("eval_im.csv",delimiter=';')

maxre = max(re)
minre = min(re)
maxim = max(im)
minim = min(im)

for i in range(len(re)):
	plt.scatter(re[i],im[i])


inp = open('input.txt')
a = inp.readline()
a = inp.readline()
a = inp.readline()
g = float(inp.readline()[2:-1])
a = float(inp.readline()[2:-1])


plt.title("real(%.4g,%.4g)  imag(%.4g,%.4g)" % (minre, maxre,minim,maxim))
if save:
	plt.savefig(name)
else:
	plt.show()



