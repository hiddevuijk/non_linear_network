import numpy as np
import matplotlib.pyplot as plt
from sys import exit

deltaEN0 = np.genfromtxt('deltaEN0.csv',delimiter=';')
deltaEN = np.genfromtxt('deltaEN.csv',delimiter=';')

plt.plot(deltaEN,color='blue',label='pop mean subtr')
plt.plot(deltaEN0,color='red',label='time avg subtr')
plt.legend()
plt.show()




