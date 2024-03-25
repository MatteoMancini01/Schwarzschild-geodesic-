# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 21:37:40 2024

@author: mmanc
"""
# =============================================================================
# Effective potential plot
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
c = 1
rs = 1
J_tilda = 2
r = np.linspace(0, 8)
Veff= (c**2 + (J_tilda/r)**2)*(1-rs/r)

discr = J_tilda**4 - 3*rs**2*c**2*J_tilda**2

r_stable  = (J_tilda**2 + np.sqrt(discr))/(rs*c**2)
r_unstable =  (J_tilda**2 - np.sqrt(discr))/(rs*c**2)

V_stable = (c**2 + (J_tilda/r_stable)**2)*(1-rs/r_stable)
V_unstable = (c**2 + (J_tilda/r_unstable)**2)*(1-rs/r_unstable)

print(r_stable, V_stable)
print(r_unstable, V_unstable)

plt.axhline(y=0.9299, color='black', label = r' $\tilde{E}^2 \approx 0.93$', linestyle = 'dashed')

plt.scatter(4.673*rs, 0.9299, color = 'purple', label = 'Perihelion shift')
plt.scatter(8.102362416407951*rs, 0.9299, color = 'purple')

plt.plot(r, Veff, color = 'black')
plt.scatter(r_stable, V_stable, color = 'blue', label = 'Stable circular orbit')
plt.scatter(r_unstable, V_unstable, color = 'red', label = 'Ustable circular orbit')
plt.title(r'Effective potential for particle with mass, $R_s = 1$')
plt.ylabel(r'$\tilde{V}^2_{eff}$')
plt.ylim(0.92,1.005)
plt.legend(loc='upper right', bbox_to_anchor=(1.45, 1))
plt.xlabel('r')

#%%
# =============================================================================
# Perihelian/Elliptical orbits
# =============================================================================
# We first need to solve the following ode using numerical integration. 
# =============================================================================
# Using RK4 Method
# =============================================================================

c = 1 # set speed of light in a vacuum equal to one
rs = 1 # schwarzchild radius 
J = 2 # angular momentum
E = 0.92999 # total energy #set to 0.92999

#starting values
phi0 = 0 # initial value for tau
r0 = 4.673*rs # initial value for r. set it to r0 = 4.673*rs


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
phi_end  = 100 # end of t values 

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
    
# here is the end of the loop we can now plot the result against the analytic
phi = np.linspace(0, 20.0, 50) # this will allow us to plot the analytic solution

plt.plot(phi_values, r_value, label = 'RK-4', color = 'blue' )

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.title('Numerical')
plt.ylabel('r')
plt.xlabel(r'$\phi$')

#%%
import numpy as np
import math
import matplotlib.pyplot as plt

# =============================================================================
# Perihelion shift, plot
# =============================================================================

r_values = [x for x in r_value if not math.isnan(x)]


r_inverse = r_values[::-1]
r_reverse = r_inverse[1:-1]

rs1 = r_values + r_reverse + r_values + r_reverse + r_values
rs2 = r_reverse + r_values + r_reverse + r_values + r_reverse

# for animation animation
rs0 = rs1 + rs2 

# shot 1
#rs0 = r_values + r_reverse + r_values

# shot 2
#rs0 =  r_values + r_reverse + r_values + r_reverse

# shot 3
#rs0 = rs1

# last shot
#rs0 = rs1 + r_reverse

a = phi_values[-1]
a_list = []
iterations = abs(len(rs0)-len(phi_values))

for _ in range(iterations):
    
    a += 0.01
    
    a_list.append(a)

phis = phi_values + a_list


if len(rs0) == len(phis):
    phi = phis
else:
    m = abs(len(rs0) - len(phis))
    phi = phis[:-m]
    

xs = []
ys = []

for i in range(len(rs0)):
    x = rs0[i]*np.cos(phi[i])
    y = rs0[i]*np.sin(phi[i])
    xs.append(x)
    ys.append(y)


# =============================================================================
# Plot of trajectory
# =============================================================================
# for a good plot set rs0 = rs1 + r_reverse

rs = 1
circle = plt.Circle((0, 0), rs, color='black', alpha=0.5, label = r'Singularity of radius $R_s$')
# Plot the trajectory in the xy-plane
plt.plot(xs, ys, label = 'perihelion orbit', color = 'blue')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.gcf().gca().add_artist(circle)

plt.xlim(-10,10)
plt.ylim(-10,10)

# Adding a text box with bullet points
textbox_content = "Paramete values:\n"
textbox_content += "• c = 1 (speed of light)\n"
textbox_content += "• $\~{J} = 2$ (angular momentum)\n"
textbox_content += "• $\~{E} = 0.92999$ (total energy)\n"
textbox_content += "• $R_s = 1$ (Schwarzschild radius)\n"
textbox_content += "• $r = 4.673 R_s$ (radial distance initial value)\n"
textbox_content += "• $\phi = 0$ (azimuthal angle initial value)\n"

plt.text(12, -10, textbox_content, bbox=dict(facecolor='white', alpha=0.5))

plt.scatter(xs[-1], ys[-1], label = "Particle's corrent position", color = 'black', s = 15 )
plt.scatter(xs[0], ys[0], label = "Particle's initial position", color = 'blue', s = 15)

plt.ylabel(r'$\frac{y}{R_s}$')
plt.xlabel(r'$\frac{x}{R_s}$')
plt.title('Perihelion shift in the xy-plane')
plt.legend(loc='upper right', bbox_to_anchor=(1.8, 1))
plt.show()
#%%

# =============================================================================
# plot of r values against phi values
# =============================================================================

plt.plot(phi, rs0)
plt.title('Numerical solution')
plt.ylabel('r')
plt.xlabel(r'$\phi$')
#plt.xlim(-1,15)

#%%

xs1 = []
ys1 = []

for i in range(len(r_values)):
    x = r_values[i]*np.cos(phi_values[i])
    y = r_values[i]*np.sin(phi_values[i])
    xs1.append(x)
    ys1.append(y)

rs = 1

circle = plt.Circle((0, 0), rs, color='black', alpha=0.5, label = r'Singularity of radius $R_s$')
# Plot the trajectory in the xy-plane
plt.plot(xs1, ys1, label = 'perihelion orbit')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.gcf().gca().add_artist(circle)

plt.xlim(-10,10)
plt.ylim(-10,10)

plt.ylabel(r'$\frac{y}{R_s}$')
plt.xlabel(r'$\frac{x}{R_s}$')
plt.title('Perihelion trajectory in the xy-plane')
plt.legend(loc='upper right', bbox_to_anchor=(1.8, 1))
plt.show()
#%%

# =============================================================================
# With animation
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib
from matplotlib import animation
# for a good animation set rs0 = rs1 + rs2

# Define the figure and axis
fig, ax = plt.subplots()

# Define the circle
rs = 1
circle = plt.Circle((0, 0), rs, color='black', alpha=0.5, label=r'Singularity of radius $R_s$')
ax.add_artist(circle)

# Plot the trajectory in the xy-plane

line, = ax.plot([], [], label='perihelion orbit', color = 'black')
ax.axhline(0, color='black', linewidth=0.5)
ax.axvline(0, color='black', linewidth=0.5)
ax.set_aspect('equal', adjustable='box')
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_ylabel(r'$\frac{y}{R_s}$')
ax.set_xlabel(r'$\frac{x}{R_s}$')
ax.set_title('Perihelion trajectory in the xy-plane')
#ax.scatter(0, 0, label = 'Particle', color = 'red', s=0)

particle, = ax.plot([], [], 'ro', label = 'Particle')  # red dot for the particle

# Function to initialize the plot
def init():
    line.set_data([], [])
    return line,

# Function to animate the plot
def animate(frame):
    x_data = xs[:frame]
    y_data = ys[:frame]
    line.set_data(x_data, y_data)
    # Update the position of the particle
    particle_x = xs[frame - 1] if frame > 0 else np.nan  # previous x coordinate
    particle_y = ys[frame - 1] if frame > 0 else np.nan  # previous y coordinate
    particle.set_data(particle_x, particle_y)
    return line, particle, 

# Create the animation
ani = FuncAnimation(fig, animate, init_func=init, frames=len(xs), interval=0, blit=True)

ax.legend(loc='upper right', bbox_to_anchor=(1.8, 1))

plt.show()
#%%

