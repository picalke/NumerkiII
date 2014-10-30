# coding: utf-8
import numerical1 as num
from defns import *
def anal(t, g):
    x=np.exp(-g/2. t)*np.cos(np.sqrt(4-g**2)/2. t)
    p=np.exp(-g/2. t)*np.sin(np.sqrt(4-g**2)/2. t)
    return x,p
def anal(t, g):
    x=np.exp(-g/2.*t)*np.cos(np.sqrt(4-g**2)/2.*t)
    p=np.exp(-g/2.*t)*np.sin(np.sqrt(4-g**2)/2.*t)
    return x,p
get_ipython().set_next_input(u't=np.arange');get_ipython().magic(u'pinfo np.arange')
t=np.arange(0., 100., 0.5)
x,p=anl(t, g) 
x,p=anal(t, g)
x
import pyplot as plt
import matplotlib.pyplot as plt
plt.figure()
plt.plot(t, x, t, p)
plt.show()
get_ipython().magic(u'save "commands.py"')
