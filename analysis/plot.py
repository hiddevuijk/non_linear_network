import numpy as np
import matplotlib.pyplot as plt

a = np.genfromtxt("deltaEN.csv",delimiter=';')
plt.plot(a)
plt.show()

