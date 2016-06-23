import numpy as np
import matplotlib.pyplot as plt
import math

def normalize(x):
	for i in reversed(range(len(x))):
		x[i] /= x[0]

def acorrArray(xt,norm=False):
	timax = xt.shape[1]
	# zero padding
	ti2 = 2*int(2**math.ceil(math.log(timax,2)))
	temp = np.lib.pad(xt,((0,0),(0,ti2-timax)),'constant',constant_values=0)

	# fourier transform
	temp = np.fft.rfft(temp,axis=1)
	#multipliy elemwise by compl.conj.
	temp = temp*np.conjugate(temp)
	#fourier back to time domain
	temp = np.fft.irfft(temp,axis=1)[:temp.shape[0],:timax]

	for i in range(temp.shape[0]):
		# devide by tmax-ti
		for j in range(temp.shape[1]):
			temp[i,j] /= temp.shape[1]-j	
		if norm == True:
			for j in reversed(range(temp.shape[1])):
				temp[i,j] /= temp[i,0]
	return temp

def acorr(x,norm=False):
	timax = x.shape[0]
	#zero padding
	ti2 = 2*int(2**math.ceil(math.log(timax,2)))
	temp = np.lib.pad(x,(0,ti2-timax),'constant',constant_values=0)
	temp = np.fft.rfft(temp)
	temp = temp*np.conjugate(temp)
	temp = np.fft.irfft(temp)[:timax]
	for j in range(temp.shape[0]):
		temp[j] /= temp.shape[0]-j
	if norm == True:
		for j in reversed(range(temp.shape[0])):
			temp[j] /= temp[0]
	return temp





