import numpy as np


def logNorm(mean,std,N):
	temp = 1+(std/mean)**2
	mu = np.log(mean/np.sqrt(temp))
	sigma = np.sqrt(np.log(temp))
	M = np.zeros((N,N))
	for i in range(N):
		for j in range(N):
			M[i,j] = np.random.lognormal(mu,sigma)
			if np.random.random_sample() > 0.5:
				M[i,j] = -1.*M[i,j]
	return M

def exp(std,N):
	M = np.random.exponential(std,(N,N))
	for i in range(N):
		for j in range(N):
			if np.random.random_sample() > 0.5:
				M[i,j] = -1.*M[i,j]
	return M


