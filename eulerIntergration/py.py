import numpy as np
import matplotlib.pyplot as plt

x = np.genfromtxt('xt.csv',delimiter=';')
plt.plot(x)
plt.show()

