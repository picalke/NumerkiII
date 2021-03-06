#consistency check -- are all variables the same type?


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
	
	
#end of consistency check

# solver. first, we define dictionary of methods

methods={'expl_midpoint': ,'impl_midpoint': , 'expl_RK4': , 'expl_euler':, 'impl_euler': }


#solver takes initial data and initial time, RHS function, dt, and max time
def solver_dt(t0, y0, f, dt, tk, method='expl_RK4', tuning='', verb=0, outform='arr'):
	
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
		tmp1, tmp2 = step(soln[i][0], soln[i][1], dt, tuning=tuning)
		soln.append([tmp1, tmp2])
	
	if v>=1: print "solved. Took "+str(nos)+" steps to solve on interval ["+str(t0)+";"+str(tk)+"]"
	if v==2: print "solve method was "+method+". Additional tweaks for this solver were:"+tuning
	if v==2: print "output form was "+outform+"(arr- numpy array such that arr[i] = [t, y1, ..., yn], list - nested list list[i]=[t, array(y)])"
	
	if outform=='arr':
		return np.array([[i[0]]+[j for j in i[1]] for i in soln])
	elif outform=='list'
		return soln
		
#solver takes initial data and initial time, RHS function, dt, and max time
def solver_nos(t0, y0, f, nos, tk, method='expl_RK4', tuning='', verb=0, outform='arr'):
	
	step=methods[method]
	#we calculate number of steps
	dt = (t0-tk)/float(nos)
		
	if v==2: print "dt "+str(dt) 
	
	#in solver_nos above code will be changed!
	
	
	soln = [[t0, y0]]
	
	for i in range(nos-1):
		tmp1, tmp2 = step(soln[i][0], soln[i][1], dt, tuning=tuning)
		soln.append([tmp1, tmp2])
	
	if v>=1: print "solved. Took "+str(nos)+" steps to solve on interval ["+str(t0)+";"+str(tk)+"]"
	if v==2: print "solve method was "+method+". Additional tweaks for this solver were:"+tuning
	if v==2: print "output form was "+outform+"(arr- numpy array such that arr[i] = [t, y1, ..., yn], list - nested list list[i]=[t, array(y)])"
	
	if outform=='arr':
		return np.array([[i[0]]+ ([a for a in i[1]] if type(i[1])==np.ndarray else [i[1]])  for i in soln])
	elif outform=='list'
		return soln
		
		
def DumpedSoln(t, g):
	A2=1.
	A1=2./np.sqrt(g**2 -4)
	x=np.exp(-g/2.*t)*(A1*np.sin(np.sqrt(g**2 -4)/2.*t)+A2*np.cos(np.sqrt(g**2 -4)/2.*t))
	p=0.5*np.exp(-g/2.*t)(np.cos(np.sqrt(g**2 -4)/2.*t)*(2-g)+np.sin(np.sqrt(g**2 -4)/2.*t)*(A1*g-2.*A2*np.sqrt(g**2 -4)/2.))
	return x, p



