# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 23:43:50 2024

@author: mmanc
"""

# =============================================================================
# Gravitational deflection
# =============================================================================
import matplotlib.pyplot as plt
import numpy as np
import math

# =============================================================================
# effective potential for photons
# =============================================================================

rs = 1
J = 2

def V(r):
    return (c**2 + (J/r)**2)*(1-rs/r)
#%%
V(3.5)

#%%

c = 1 # set speed of light in a vacuum equal to one
rs = 1 # schwarzchild radius 
J = 2 # angular momentum
E = V(3.5) + 0.01 # total energy #set to V(3) + 0.01

#starting values
phi0 = 0 # initial value for tau
r0 = 3.159*rs # initial value for r. set it to r0 = 2.93397443*rs

# now we are going to define f(y,t) our ODE

def f(phi, r):
    return (r**2/J) * np.sqrt(E - (1 - rs/r)*(c**2 + (J/r)**2))# remember this is plus or minus

# Now we are going to set few lists:

phi_values = [] # will collect all the values of tau
r_value = [] # will collect all the values of r at proper time tau

phi_old = phi0
r_old = r0 # used for iterating procedure 

phi_values.append(phi_old)
r_value.append(r_old) # appending old values to lists

# let delta in the RK-4 method be rappresented by h, the smaller
# we chose h to be the lower will the error be (h is the step size):
h = 0.01 # chose h = 0.01 for accurate result
phi_end  = 20*np.pi # end of t values 

Nsteps = np.int(np.floor((phi_end-phi0)/h)) # this will allow us to define a range in 
# our for loop
# now we can start our loop
for i in range(0, Nsteps):
    
    phi_now = phi_old + h
    
    k1 = f(phi_old,r_old)
    k2 = f(phi_old + 0.5*h, r_old + 0.5*k1*h)
    k3 = f(phi_old + 0.5*h, r_old + 0.5*k2*h)
    k4 = f(phi_old + h, r_old + k3*h)
    
    r_now = r_old + h*(k1 + 2*k2 + 2*k3 + k4)/6
    
    phi_values.append(phi_now)
    r_value.append(r_now)
    
    phi_old = phi_now
    r_old = r_now
    

plt.plot(phi_values, r_value, label = 'RK-4', color = 'black' )

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.title('Numerical solution')
plt.ylabel('r')
plt.xlabel(r'$\phi$')

#%%
phi_negative = []

for i in phi_values:
    phi_negative.append(-i)
phi_value = phi_negative + phi_values

rvalue = r_value + r_value

plt.plot(phi_value, rvalue, color = 'black', label = 'RK4')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.title('Numerical solution')
plt.ylabel('r')
plt.xlabel(r'$\phi$')


#%%
# convering to cartesian
# convering to cartesian
x = rvalue*np.cos(phi_value)
y = rvalue*np.sin(phi_value)

circle = plt.Circle((0, 0), rs, color='black', alpha=0.5, label = r'Singularity of radius $R_s$')

# Plot the trajectory in the xy-plane
plt.plot(x, y, label = 'Gravitational deflection', color = 'blue')


plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.gcf().gca().add_artist(circle)


# Adding a text box with bullet points
textbox_content = "Paramete values:\n"
textbox_content += "• c = 1 (speed of light)\n"
textbox_content += "• $\~{J} = 2$ (angular momentum)\n"
textbox_content += "• $\~{E}^2 \u2248 0.95 $ (total energy)\n"
textbox_content += "• $R_s = 1$ (Schwarzschild radius)\n"
textbox_content += "• $r \u2248 3.2 R_s$ (radial distance initial value)\n"
textbox_content += "• $\phi = 0$ (azimuthal angle initial value)\n"

plt.text(22, -12, textbox_content, bbox=dict(facecolor='white', alpha=0.5))

d = 20

plt.xlim(-d,d)
plt.ylim(-d,d)

plt.ylabel(r'$\frac{y}{R_s}$')
plt.xlabel(r'$\frac{x}{R_s}$')
plt.title('Gravitational deflection of mass particles')
plt.legend(loc='upper right', bbox_to_anchor=(1.8, 1))
plt.show()

