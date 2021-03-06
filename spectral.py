#! /usr/bin/python

import numpy as np
import scipy as sp
import scipy.optimize as opt
import numerical1 as num
from datetime import datetime


#trick for replacing U-k with Uk* 
def inv(x):
	def d(i):
	        if i>=0 and i<len(x):
	            return x[i]
	        elif i<0 and -i<len(x):
	            return np.conj(x[-i])
	        else:
	            return 0.+0.j
	return d

def E(y):
	uc=inv(y[:len(y)/2])
	r=range(-len(y[:len(y)/2])+1, len(y[:len(y)/2]))
	enon=np.sum([[[uc(-k)*uc(-l)*uc(m)*uc(k+l-m) for m in r] for l in r] for k in r])
	elin=np.pi*np.sum([np.abs(y[len(y/2):][k]) + k**2*np.abs(y[:len(y/2)][k]) for k in range(len(y)/2)])
	return enon+elin

def TimeDeriv(t, y):
	uc=inv(y[:len(y)/2])
	le=len(y[:len(y)/2])
	pdot=-np.array([i**2*uc(i) for i in range(le)]) - np.array([np.sum([[uc(l)*uc(k)*uc(n-k-l) for l in range(-le+1, le)] for k in range(-le+1, le)]) for n in range(le)])	
	return np.reshape([y[len(y)/2:], pdot], 2*le)

#def TimeDeriv(y):
#	
#	return

#spectral solver takes initial data in fourier
def spectralSolve(u0, ut0, tmax, ntime, nsave ,method='expl_RK4', filename="./spectral.dat", file2="./spectralEnergy.dat"):
	dt = tmax/float(ntime)
	t=0.
	y=np.reshape([u0, ut0], len(u0)*2)
	print dt
	plik=open(filename, "w")
	plik2=open(file2, "w")
	#add some header to the file at some point plax
	plik2.write(str(t)+'\t'+str(E(y))+'\n')
	ynew=np.zeros(len(y), dtype=np.complex)
	string = str(t)+'\t'+''.join([str(el)+"\t" for el in y[:len(ynew)/2]])+'\n'
	plik.write(string)
	for st in range(ntime):
		tstart=datetime.now()
		t, ynew = num.methods[method](TimeDeriv, y, t, dt)
		tend=datetime.now()
		if st%nsave==0:
			string = str(t)+'\t'+''.join([str(el)+"\t" for el in y[:len(y)/2]])+'\n'
			plik.write(string)
			plik2.write(str(t)+'\t'+str(E(y))+'\n')
		ynew, y = y, ynew
	
	string = str(t)+'\t'+''.join([str(el)+"\t" for el in y[:len(ynew)/2]])+'\n'
	plik.write(string)
	plik2.write(str(t)+'\t'+str(E(y))+'\n')
	plik.close()
	plik2.close()
	return 0.
	

	
	
	
