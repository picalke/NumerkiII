import numpy as np
import numerical1 as num

g=0.000001

def osc(t,y):
	global g	
	return np.dot(np.array([[-1*g, -1],[1,0]]), y)

y0=np.array([0,1])

def DumpedSoln(t, g):
	A2=1.
	A1=2./np.sqrt(-g**2+4)
	x=np.exp(-1*g/2.*t)*(A1*np.sin(np.sqrt(-g**2+4)/2.*t)+A2*np.cos(np.sqrt(-g**2 +4)/2.*t))
	p=0.5*np.exp(-1*g/2.*t)*(np.cos(np.sqrt(-g**2+4)/2.*t)*(2-g)+np.sin(np.sqrt(-g**2+4)/2.*t)*(A1*g-2.*A2*np.sqrt(-g**2+4)/2.))
	return x, p


def anal(t, g):
    x=np.exp(-g/2.*t)*np.cos(np.sqrt(4-g**2)/2.*t)
    p=np.exp(-g/2.*t)*np.sin(np.sqrt(4-g**2)/2.*t)
    return x,p
    
import matplotlib.pyplot as plt

sol=num.solver_dt(0., y0, osc, 0.05, 100, method='expl_euler')
sol2=num.solver_dt(0., y0, osc, 0.05, 100, method='impl_midpoint')



