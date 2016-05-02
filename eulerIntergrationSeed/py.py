import numpy as np
import matplotlib.pyplot as plt

x = np.genfromtxt('deltaN_123456789.csv',delimiter=';')
plt.plot(x)
plt.show()

