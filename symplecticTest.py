import numpy as np
import numerical1 as num
import matplotlib.pyplot as plt

def gradU(q):
	return q
	
def gradT(p):
	return p

y0=np.array([0,1])

def H(p,q):
	return p**2+q**2

def osc(t):
	p=np.sin(t)
	q=np.cos(t)
	return p, q

sol = num.solver_nos(0., y0, (gradU,gradT), 2000, 20*np.pi, method='neri')
