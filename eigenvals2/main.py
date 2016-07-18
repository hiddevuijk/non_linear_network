import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
from dist import *

parser = argparse.ArgumentParser()

parser.add_argument('--size','-N',dest='N',type=int,default=500)
parser.add_argument('--mu','-m',dest='mu',type=float,default=0.0)
parser.add_argument('--sigma','-s',dest='s',type=float,default=1.)
parser.add_argument('--dist','-d',dest='dist',type=str,default='normal')
parser.add_argument('--seed',dest='seed',type=int,default=0)
parser.add_argument('--rmMax',dest='rmMax',action='store_true',default=False)


args = parser.parse_args()
rmMax  = args.rmMax
N = args.N
mu = args.mu
s = args.s
dist = args.dist
seed = args.seed

distOpts =  ['normal','lognormal','exp']

if dist not in distOpts:
	print dist, " is not a valid option for dist."
	print "Valid options are: "
	for d in distOpts:
		print '\t',d
	sys.exit(1)

# create random matrix M
#if seed is not 0, set random generatro seed to seed
if seed != 0:
	np.random.seed(seed)

if dist == 'normal':
	M = np.random.normal(mu,s,(N,N))
	
elif dist == 'lognormal':
	M = logNorm(mu,s,N)
elif dist == 'exp':
	M = exp(s,N) 

# get eigenvalues
eVal, eVec = np.linalg.eig(M)

if rmMax == True:
	eVal = np.sort(eVal)[:-1]


# make plot

ax1 = plt.subplot2grid((6,4),(0,0),rowspan=4,colspan=4)
ax2 = plt.subplot2grid((6,4),(4,0),rowspan=2,colspan=4)

# plot real,imag eVal

ax1.scatter(eVal.real,eVal.imag)
ax1.axhline(max(eVal.imag),color='black',linestyle='--')
ax1.axhline(min(eVal.imag),color='black',linestyle='--')
ax1.axvline(max(eVal.real),color='black',linestyle='--')
ax1.axvline(min(eVal.real),color='black',linestyle='--')
ax1.set_aspect('equal')
title = 'std = '+str(M.std())
ax1.set_title(title)

xmax = max([max(eVal.real),max(eVal.imag)])
ax2.hist(abs(eVal),normed=True)
ax2.set_xlim([0,xmax])

plt.tight_layout()
plt.show()





