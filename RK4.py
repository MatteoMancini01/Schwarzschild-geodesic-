# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

# the given ODE is y' = -y, we want to solve this numerically using
#RK-4 method, to verify that the method is accurate we can do a plot,
# of the analytic solution and the numerical solution
# Analytic => y(t) = A*exp(-t) general sol, A = const

#defining analytic sol
A = 1
#starting values
t0 = 0
y0 = 1

def Y(t):
    return A*np.exp(-t)

# above we have the analytic solution fix A = 1, so that y(0) = 1.

# now we are going to define f(y,t) our ODE

def f(t,y):
    return -y

# Now we are going to set few lists:

ts = [] # will collect all the values of t
ys = [] # will collect all the values of y at time t

told = t0
yold = y0 # used for iterating procedure 

ts.append(told)
ys.append(yold) # appending old values to lists

# let delta in the RK-4 method be rappresented by h, the smaller
# we chose h to be the lower will the error be (h is the step size):
h = 0.01 # chose h = 0.5 for now
tend  = 10 # end of t values 

Nsteps = np.int(np.floor((tend-t0)/h)) # this will allow us to define a range in 
# our for loop
# now we can start our loop
for i in range(0, Nsteps):
    
    tnow = told + h
    
    k1 = f(told,yold)
    k2 = f(told + 0.5*h, yold + 0.5*k1*h)
    k3 = f(told + 0.5*h, yold + 0.5*k2*h)
    k4 = f(told + h, yold + k3*h)
    
    ynow = yold + h*(k1 + 2*k2 + 2*k3 + k4)/6
    
    ts.append(tnow)
    ys.append(ynow)
    
    told = tnow
    yold = ynow
    
# here is the end of the loop we can now plot the result against the analytic
t = np.linspace(0, 10.0, 50) # this will allow us to plot the analytic solution
analytic = Y(t)
plt.plot(ts, ys,label = 'RK-4', color = 'blue' )
plt.plot(t,analytic, label = 'analytic', color = 'red')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.title('Numerical vs analytic')
plt.ylabel('y')
plt.xlabel('t')
plt.grid()
# By testing the output of the code using different steps we see that the accuracy,
# of the RK4 increases when h is small
#%%
# =============================================================================
# Now compare the above with the plot obtained using scipy.integrate
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# define function that returns dy/dt

def model(y,t):
    dydt = -y
    return dydt

#initial coordinate
y0 = 1
#time points
t = np.linspace(0,10)

# solve ODE model
y = odeint(model, y0, t)

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()

# as we can see the two plots are the same. To solve ODEs that model motion of a paticle
# We are going to use this python package.






