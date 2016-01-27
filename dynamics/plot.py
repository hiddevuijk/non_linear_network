import matplotlib.pyplot as plt
import numpy as np


xt = np.genfromtxt("x.csv",delimiter=';')
t = np.genfromtxt("t.csv",delimiter=';')

N= xt.shape[0]
dt=xt.shape[1]

N = min(N,3)

for i in range(N):
	plt.plot(t,xt[i])

plt.show()





