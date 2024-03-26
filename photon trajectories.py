# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 21:47:31 2024

@author: mmanc
"""

import matplotlib.pyplot as plt
import numpy as np
import math

# =============================================================================
# effective potential for photons
# =============================================================================

rs = 1
J = 2

def V(r):
    return (1 - rs/r)*(J/r)**2

r = np.linspace(0, 10)

plt.axhline(y=V(1.5), color='black', label = r'Total energy  $E^2_0 \approx 0.6$', linestyle = 'dashed' )
plt.axhline(y=V(3)+0.01, color='black', label = r'Total energy  $E^2_0 \approx 0.4$', linestyle = 'dashdot' )
plt.plot(r, V(r), color = 'black')
plt.scatter(1.5, V(1.5), label = "circular orbit 'unstable'")
plt.scatter(2.93, V(3)+0.01, label = "light deflection", color = 'red')
plt.title(r'Effective potential of photons')
plt.xlabel('r')
plt.ylabel(r'$V_{eff}$(r)')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.ylim(0,0.7)
plt.show()
#%%
# =============================================================================
# RK4 Method
# =============================================================================

c = 1 # set speed of light in a vacuum equal to one
rs = 1 # schwarzschild radius 
J = 2 # angular momentum
E = V(1.5) # total energy #set to 0.92999

#starting values
phi0 = 0 # initial value for tau
r0 = 1.5*rs # initial value for r. set it to r0 = 4.673*rs


# now we are going to define f(y,t) our ODE

def f(phi, r):
    return np.sqrt(E*(r**2/J)**2 - r**2 + rs*r)# remember this is plus or minus

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
phi_end  = 2*np.pi # end of t values 

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
    

plt.plot(phi_values, r_value, label = 'RK-4', color = 'blue' )

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.title('Numerical solution')
plt.ylabel('r')
plt.xlabel(r'$\phi$')

#%%

# =============================================================================
# 
# =============================================================================

xs1 = []
ys1 = []

for i in range(len(r_value)):
    x = r_value[i]*np.cos(phi_values[i])
    y = r_value[i]*np.sin(phi_values[i])
    xs1.append(x)
    ys1.append(y)


rs = 1

circle = plt.Circle((0, 0), rs, color='black', alpha=0.5, label = r'Singularity of radius $R_s$')

# Plot the trajectory in the xy-plane
plt.plot(xs1, ys1, label = 'Unstable circular orbit', color = 'orange')


plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.gcf().gca().add_artist(circle)


# Adding a text box with bullet points
textbox_content = "Paramete values:\n"
textbox_content += "• c = 1 (speed of light)\n"
textbox_content += "• $J = 2$ (angular momentum)\n"
textbox_content += "• $E^2 \u2248 0.6 $ (total energy)\n"
textbox_content += "• $R_s = 1$ (Schwarzschild radius)\n"
textbox_content += "• $r = 3/2 R_s$ (radial distance)\n"
#textbox_content += "• $\phi = 0$ (azimuthal angle initial value)\n"

plt.text(3.5, -2, textbox_content, bbox=dict(facecolor='white', alpha=0.5))

d = 3

plt.xlim(-d,d)
plt.ylim(-d,d)

plt.ylabel(r'$\frac{y}{R_s}$')
plt.xlabel(r'$\frac{x}{R_s}$')
plt.title('Unstable circular orbit photons')
plt.legend(loc='upper right', bbox_to_anchor=(1.8, 1))
plt.show()
#%%
# =============================================================================
# Gravitational deflection
# =============================================================================

c = 1 # set speed of light in a vacuum equal to one
rs = 1 # schwarzchild radius 
J = 2 # angular momentum
E = V(3) + 0.01 # total energy #set to V(3) + 0.01

#starting values
phi0 = 0 # initial value for tau
r0 = 2.93397443*rs # initial value for r. set it to r0 = 2.93397443*rs


# now we are going to define f(y,t) our ODE

def f(phi, r):
    return np.sqrt(E*(r**2/J)**2 - r**2 + rs*r)# remember this is plus or minus

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
phi_end  = 10*np.pi # end of t values 

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
    

plt.plot(phi_values, r_value, label = 'RK-4', color = 'blue' )

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.title('Numerical')
plt.ylabel('r')
plt.xlabel(r'$\phi$')


#%%

xs = []
ys = []

for i in range(len(r_value)):
    x = r_value[i]*np.cos(phi_values[i])
    y = r_value[i]*np.sin(phi_values[i])
    xs.append(x)
    ys.append(y)

xs0 = []
ys0 = []
for j in range(len(r_value)):
    x0 = r_value[j]*np.cos(phi_values[j])
    y0 = - r_value[j]*np.sin(phi_values[j])
    xs0.append(x0)
    ys0.append(y0)

xs1 = xs0 + xs
ys1 = ys0 + ys

rs = 1

circle = plt.Circle((0, 0), rs, color='black', alpha=0.5, label = r'Singularity of radius $R_s$')

# Plot the trajectory in the xy-plane
plt.plot(xs1, ys1, label = 'Gravitational deflection', color = 'orange')


plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.gcf().gca().add_artist(circle)


# Adding a text box with bullet points
textbox_content = "Paramete values:\n"
textbox_content += "• c = 1 (speed of light)\n"
textbox_content += "• $J = 2$ (angular momentum)\n"
textbox_content += "• $E^2 \u2248 0.3 $ (total energy)\n"
textbox_content += "• $R_s = 1$ (Schwarzschild radius)\n"
textbox_content += "• $r \u2248 2.93 R_s$ (radial distance initial value)\n"
textbox_content += "• $\phi = 0$ (azimuthal angle initial value)\n"

plt.text(18, -10, textbox_content, bbox=dict(facecolor='white', alpha=0.5))

d = 15

plt.xlim(-d,d)
plt.ylim(-d,d)

plt.ylabel(r'$\frac{y}{R_s}$')
plt.xlabel(r'$\frac{x}{R_s}$')
plt.title('Gravitational deflection of light')
plt.legend(loc='upper right', bbox_to_anchor=(1.8, 1))
plt.show()

















