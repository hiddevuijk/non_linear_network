import matplotlib.pyplot as plt
import numpy as np

r = 15

re = np.genfromtxt("eval_re.csv",delimiter=';')
im = np.genfromtxt("eval_im.csv",delimiter=';')


for i in range(len(re)):
	plt.scatter(re[i],im[i])
#plt.ylim([-1*r,r])
#plt.xlim([-1*r,r])
plt.show()



