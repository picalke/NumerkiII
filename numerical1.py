#! /usr/bin/python

import numpy as np
import scipy as sp
import scipy.optimize as opt

#inner method! do not use in general unless you know syntax!
def vec(a):
	return np.reshape(a, 2*len(a[0]))

#one step of method. Full integrator in another function
#f(t, y)--LHS function in standard form, y -- current value, t -- current time, dt -- step
#
def expl_euler_step(f, y, t, dt, tuning=[]):
	if type(y)!=type(f(t,y)):
		print "solver error 0: type returned by RHS not matching type of data!"
		return float('NaN'),float('NaN')
	elif type(y)==np.ndarray:
		if len(y)!=len(f(t,y)):
			print "solver error 1: length of array returned by RHS not matching size of data array!"
			return float('NaN'),float('NaN')
	
	elif type(y)!=float and type(f(t,y))!=float and type(y)!=np.float64 and type(f(t,y))!=np.float64:
		print "solver error 3: Type of data or RHS does not match desired format (float or numpy.float64)"
		return float('NaN'),float('NaN')

	return t+dt, y+dt*f(t,y)


####
def expl_mpm_step(f, y, t, dt, tuning=[]):
	if type(y)!=type(f(t,y)):
		print "solver error 0: type returned by RHS not matching type of data!"
		return float('NaN'),float('NaN')
	elif type(y)==np.ndarray:
		if len(y)!=len(f(t,y)):
			print "solver error 1: length of array returned by RHS not matching size of data array!"
			return float('NaN'),float('NaN')
	
	elif type(y)!=float and type(f(t,y))!=float and type(y)!=np.float64 and type(f(t,y))!=np.float64:
		print "solver error 3: Type of data or RHS does not match desired format (float or numpy.float64)"
		return float('NaN'),float('NaN')


	k1=f(t, y)
	k2=f(t+0.5*dt, y+0.5*dt*k1)
	yn = y+dt*k2
	return t+dt, yn
####
def impl_mpm_step(f, y, t, dt, tuning=['default']):
	if type(y)!=type(f(t,y)):
		print "solver error 0: type returned by RHS not matching type of data!"
		return float('NaN'),float('NaN')
	elif type(y)==np.ndarray:
		if len(y)!=len(f(t,y)):
			print "solver error 1: length of array returned by RHS not matching size of data array!"
			return float('NaN'),float('NaN')
	
	elif type(y)!=float and type(f(t,y))!=float and type(y)!=np.float64 and type(f(t,y))!=np.float64:
		print "solver error 3: Type of data or RHS does not match desired format (float or numpy.float64)"
		return float('NaN'),float('NaN')

	#by defauld we use Newton method with Krylow approx of jacobian
	#maybe later some if statements and options, who knows
	sol = opt.fsolve(lambda yn: y+dt*f(t+0.5*dt, 0.5*(y+yn))-yn, y)[0] if type(y)!=np.ndarray else opt.fsolve(lambda yn: y+dt*f(t+0.5*dt, 0.5*(y+yn))-yn, y)
	return t+dt, sol
	
######
def expl_rk4_step(f, y, t, dt, tuning=['default']):
	if type(y)!=type(f(t,y)):
		print "solver error 0: type returned by RHS not matching type of data!"
		return float('NaN'),float('NaN')
	elif type(y)==np.ndarray:
		if len(y)!=len(f(t,y)):
			print "solver error 1: length of array returned by RHS not matching size of data array!"
			return float('NaN'),float('NaN')
	
	elif type(y)!=float and type(f(t,y))!=float and type(y)!=np.float64 and type(f(t,y))!=np.float64:
		print "solver error 3: Type of data or RHS does not match desired format (float or numpy.float64)"
		return float('NaN'),float('NaN')

	k1=f(t, y)
	k2=f(t+0.5*dt, y+dt*0.5*k1)
	k3=f(t+0.5*dt, y+dt*0.5*k2)
	k4=f(t+dt, y+dt*k1)
	yn = y+dt*(0.16666666666666666*k1+0.3333333333333333*k2+0.3333333333333333*k3+0.16666666666666666*k4)
	return t+dt, yn


def expl_sym_euler_step(f, y, t, dt, tuning=['default']):
	if type(y)!=type(f(t,y)):
		print "solver error 0: type returned by RHS not matching type of data!"
		return float('NaN'),float('NaN')
	elif type(y)==np.ndarray:
		if len(y)!=len(f(t,y)):
			print "solver error 1: length of array returned by RHS not matching size of data array!"
			return float('NaN'),float('NaN')
		if len(y)%2!=0 or len(f(t,y))%2!=0:
			print "symplectic method only valid for Hamiltonian system: such system needs to be at least two-dimensional and of even dimension, which is not the case here!"
			return float('NaN'), float('NaN')
	else:
		print "for symplectic method you need to have a Hamiltonian system which means f returns at least 2-dimensional vector!"
		return float('NaN'),float('NaN')
	#expecting that f = J^-1 grad H(p,q)!!!!!and H = T(p)+U(q)!!!!
	halba=len(y)/2
	pn=y[:halba]+dt*f(t,y)[:halba]
	tmp_y=np.reshape([pn, y[halba:]])
	#Ty, mosz recht, wez dwie!
	qn=y[halba:]+dt*f(t,tmp_y)[halba:]
	yn = np.reshape([pn,qn], halba*2)
	return t+dt, yn
	
	 
def expl_sym_verlet_step(f, y, t, dt, tuning=['default']):
	if type(y)!=type(f(t,y)):
		print "solver error 0: type returned by RHS not matching type of data!"
		return float('NaN'),float('NaN')
	elif type(y)==np.ndarray:
		if len(y)!=len(f(t,y)):
			print "solver error 1: length of array returned by RHS not matching size of data array!"
			return float('NaN'),float('NaN')
		if len(y)%2!=0 or len(f(t,y))%2!=0:
			print "symplectic method only valid for Hamiltonian system: such system needs to be at least two-dimensional and of even dimension, which is not the case here!"
			return float('NaN'), float('NaN')
	else:
		print "for symplectic method you need to have a Hamiltonian system which means f returns at least 2-dimensional vector!"
		return float('NaN'),float('NaN')
	#expecting that f = -J^-1 grad H(p,q)!!!!!and H = T(p)+U(q)!!!!
	halba=len(y)/2
	tmp_p = y[:halba] +0.5*dt*f(t,y)[:halba]
	tmp_q = 0
	print "UNFINISHED"
	return float('NaN'), float('NaN')
	
def expl_neri_step(f, y, t, dt, tuning=['default']):
	if type(y)!=type(f(t,y)):
		print "solver error 0: type returned by RHS not matching type of data!"
		return float('NaN'),float('NaN')
	elif type(y)==np.ndarray:
		if len(y)!=len(f(t,y)):
			print "solver error 1: length of array returned by RHS not matching size of data array!"
			return float('NaN'),float('NaN')
		if len(y)%2!=0 or len(f(t,y))%2!=0:
			print "symplectic method only valid for Hamiltonian system: such system needs to be at least two-dimensional and of even dimension, which is not the case here!"
			return float('NaN'), float('NaN')
	else:
		print "for symplectic method you need to have a Hamiltonian system which means f returns at least 2-dimensional vector!"
		return float('NaN'),float('NaN')
	#expecting that f = -J^-1 grad H(p,q)!!!!!and H = T(p)+U(q)!!!!
	halba=len(y)/2
	c=np.array([0.15337808457192392,-0.03986574628422823,-0.03986574628422823,0.15337808457192392,])
	d=np.array([0.30675616914384785, -0.3864876617123043,0.30675616914384785,0.])
	tmp_new = np.array([0.]*len(y))
	tmp_old = y.copy()
	for i in range(len(c)):
		tmp_new[halba:]=tmp_old[halba:]+c[i]*dt*f(tmp_old)[halba:]
		tmp_new[:halba]=tmp_old[:halba]+d[i]*dt*f(tmp_new)[:halba]
		tmp=tmp_new[:]
		tmp_new=tmp_old[:]
		tmp_old=tmp[:]
	return t+dt, tmp_old

def expl_suzuki_step(f, y, t, dt, tuning=['default']):
	if type(y)!=type(f(t,y)):
		print "solver error 0: type returned by RHS not matching type of data!"
		return float('NaN'),float('NaN')
	elif type(y)==np.ndarray:
		if len(y)!=len(f(t,y)):
			print "solver error 1: length of array returned by RHS not matching size of data array!"
			return float('NaN'),float('NaN')
		if len(y)%2!=0 or len(f(t,y))%2!=0:
			print "symplectic method only valid for Hamiltonian system: such system needs to be at least two-dimensional and of even dimension, which is not the case here!"
			return float('NaN'), float('NaN')
	else:
		print "for symplectic method you need to have a Hamiltonian system which means f returns at least 2-dimensional vector!"
		return float('NaN'),float('NaN')
	#expecting that f = -J^-1 grad H(p,q)!!!!!and H = T(p)+U(q)!!!!
	halba=len(y)/2
	c=np.array([ 0.20724476,  0.41448951, -0.12173427, -0.12173427,  0.41448951,
        0.20724476])
	d=np.array([ 0.41448951,  0.41448951, -0.65795805,  0.41448951,  0.41448951,  0.        ])
	tmp_new = np.array([0.]*len(y))
	tmp_old = y.copy()
	for i in range(len(c)):
		tmp_new[halba:]=tmp_old[halba:]+c[i]*dt*f(tmp_old)[halba:]
		tmp_new[:halba]=tmp_old[:halba]+d[i]*dt*f(tmp_new)[:halba]
		tmp=tmp_new[:]
		tmp_new=tmp_old[:]
		tmp_old=tmp[:]
	return t+dt, tmp_old	
#any integrator

# solver. first, we define dictionary of methods

methods={'expl_midpoint':expl_mpm_step ,'impl_midpoint':impl_mpm_step , 'expl_RK4':expl_rk4_step , 'expl_euler':expl_euler_step, 'expl_symplectic_euler':expl_sym_euler_step, 'neri':expl_neri_step, 'suzuki':expl_suzuki_step}


#solver takes initial data and initial time, RHS function, dt, and max time
def solver_dt(t0, y0, f, dt, tk, method='expl_RK4', tuning='', v=0, outform='arr'):
	
	step=methods[method]
	#we calculate number of steps
	if np.abs(t0-tk/dt)-np.floor(np.abs(t0-tk/dt)) >= 0.5:
		nos = int(np.ceil(np.abs(t0-tk)/dt))
	else:
		nos = int(np.floor(np.abs(t0-tk)/dt))
		
	if v==2: print "Number of steps I`ll take is "+str(nos) 
	
	#in solver_nos above code will be changed!
	
	
	soln = [[t0, y0]]
	
	for i in range(nos-1):
		tmp1, tmp2 = step(f, soln[i][1], soln[i][0], dt, tuning=tuning)
		soln.append([tmp1, tmp2])
	
	if v>=1: print "solved. Took "+str(nos)+" steps to solve on interval ["+str(t0)+";"+str(tk)+"]"
	if v==2: print "solve method was "+method+". Additional tweaks for this solver were:"+tuning
	if v==2: print "output form was "+outform+"(arr- numpy array such that arr[i] = [t, y1, ..., yn], list - nested list list[i]=[t, array(y)])"
	
	if outform=='arr':
		return np.array([[i[0]]+ ([a for a in i[1]] if type(i[1])==np.ndarray else [i[1]])  for i in soln])
	elif outform=='list':
		return soln
		
#solver takes initial data and initial time, RHS function, dt, and max time
def solver_nos(t0, y0, f, nos, tk, method='expl_RK4', tuning='', v=0, outform='arr'):
	
	step=methods[method]
	#we calculate number of steps
	dt = (t0-tk)/float(nos)
		
	if v==2: print "dt "+str(dt) 
	
	#in solver_nos above code will be changed!
	
	
	soln = [[t0, y0]]
	
	for i in range(nos-1):
		tmp1, tmp2 = step(f, soln[i][1], soln[i][dt], 0, tuning=tuning)
		soln.append([tmp1, tmp2])
	
	if v>=1: print "solved. Took "+str(nos)+" steps to solve on interval ["+str(t0)+";"+str(tk)+"]"
	if v==2: print "solve method was "+method+". Additional tweaks for this solver were:"+tuning
	if v==2: print "output form was "+outform+"(arr- numpy array such that arr[i] = [t, y1, ..., yn], list - nested list list[i]=[t, array(y)])"
	
	if outform=='arr':
		return np.array([[i[0]]+ ([a for a in i[1]] if type(i[1])==np.ndarray else [i[1]])  for i in soln])
	elif outform=='list':
		return soln
		




	
