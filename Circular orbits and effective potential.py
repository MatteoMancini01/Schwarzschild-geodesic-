# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 21:41:10 2024

@author: mmanc
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
# =============================================================================
# Effective potential for massive objects
# =============================================================================
c = 1 # set the speed of light equal to 1
rs = 1.3# also set the Schwarzchild radius to 1
#angular momentum
L_tilda0 = 1
L_tilda1 = 2
L_tilda2 = 8
L_tilda3 = 10

r = np.linspace(0,10)

Vetilda0 = (c**2 + (L_tilda0/r)**2)*(1-rs/r)
Vetilda1 = (c**2 + (L_tilda1/r)**2)*(1-rs/r)
Vetilda2 = (c**2 + (L_tilda2/r)**2)*(1-rs/r)
Vetilda3 = (c**2 + (L_tilda3/r)**2)*(1-rs/r)

plt.plot(r, Vetilda0, color = 'red', label = r'$\tilde{J} = 1$')
plt.plot(r, Vetilda1, color = 'orange', label = r'$\tilde{J} = 2$')
plt.plot(r, Vetilda2, color = 'blue', label = r'$\tilde{J} = 8$')
plt.plot(r, Vetilda3, color = 'purple', label = r'$\tilde{J} = 10$')
plt.title(r'Effective potential for particle with mass for $R_s = 1.3$')
plt.ylabel(r'$\tilde{V}^2_{eff}$')
plt.ylim(-5,10)
plt.xlabel('r')
plt.legend()

#%%
# =============================================================================
# Effective potential for circular orbit (stable/unstable)
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




plt.plot(r, Veff, color = 'black')
plt.scatter(r_stable, V_stable, color = 'blue', label = 'Stable circular orbit')
plt.axhline(y=V_stable, color='black', label = r' $\tilde{E}_2^2 \approx 0.926$', linestyle = 'dashdot')

plt.scatter(r_unstable, V_unstable, color = 'red', label = 'Ustable circular orbit')
plt.axhline(y=V_unstable, color='black', label = r' $\tilde{E}_1^2 \approx 1.0$', linestyle = 'dashed')
plt.title(r'Effective potential for particle with mass, $R_s = 1$')
plt.ylabel(r'$\tilde{V}^2_{eff}$')
plt.ylim(0.92,1.005)
plt.legend(loc='upper right', bbox_to_anchor=(1.45, 1))
plt.xlabel('r')

#%%
# =============================================================================
# Circular orbits stable and unstable
# =============================================================================
# If we set r to be constant, this leads to the simplest case of orbit; the circular orbit
# one can find the value of r for a circular orbit by taking the derivative with respect to r,
# of the effective potential leading to the following:

J_tilda = 2 # set it to 2
rs = 1 # set it to 1
c = 1

discr = J_tilda**4 - 3*rs**2*c**2*J_tilda**2

circle = plt.Circle((0, 0), rs, color='black', alpha=0.5, label = r'Singularity of radius $R_s$')

singcircle = plt.Circle((0, 0), rs, color='black', alpha=0.5, label = r'Singularity of radius $R_s$')

if discr >=0:#when discriminatnt is positive
    
    r_stable  = (J_tilda**2 + np.sqrt(discr))/(rs*c**2)

    r_unstable =  (J_tilda**2 - np.sqrt(discr))/(rs*c**2)

    period = np.linspace(0, 2*np.pi)


    x_value_stable = r_stable*np.cos(period)
    y_value_stable = r_stable*np.sin(period)

    x_value_unstable = r_unstable*np.cos(period)
    y_value_unstable = r_unstable*np.sin(period)
    
    plt.plot(x_value_unstable, y_value_unstable, color = 'red', label = 'Unstable circular orbit')
    plt.scatter(r_unstable*np.cos(3*np.pi/4), r_unstable*np.sin(3*np.pi/4), color = 'red', label = 'Particle1')

    plt.plot(x_value_stable, y_value_stable, color = 'blue', label = 'Stable circular orbit')
    plt.scatter(r_stable*np.cos(np.pi/4), r_stable*np.sin(np.pi/4), color = 'blue', label = 'Particle2')

    plt.gcf().gca().add_artist(circle)

else:# when discriminant is negative
    
    plt.gcf().gca().add_artist(singcircle)
    plt.scatter((rs-rs/4)*np.cos(np.pi/4), (rs-rs/4)*np.sin(np.pi/4), color = 'red', label = 'Particle1')
    plt.scatter((rs-rs/3)*np.cos(3*np.pi/4), (rs-rs/3)*np.sin(3*np.pi/4), color = 'blue', label = 'Particle2')
    plt.ylim(-rs-5, rs+5)
    plt.xlim(-rs-5, rs+5)

#plt.text(-8, -9.5, r'Here the Schwarzchild radius is set as $R_s = 1$ and angular momentum is $\tilde{J} = 2$')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.gca().set_aspect('equal', adjustable='box')

plt.title('Circular orbits for the Schwarzschild geodesic')
plt.ylabel(r'$\frac{y}{R_s}$')
plt.xlabel(r'$\frac{x}{R_s}$')
plt.legend(loc='upper right', bbox_to_anchor=(1.8, 1))

# Adding a text box with bullet points
textbox_content = "Paramete values:\n"
textbox_content += "• c = 1 (speed of light)\n"
textbox_content += "• $R_s = 1$ (Schwarzschild radius)\n"
textbox_content += "• $\~{J} = 2$ (angular momentum)\n"

plt.text(8, -5, textbox_content, bbox=dict(facecolor='white', alpha=0.5))

plt.show()

print(r_stable, r_unstable)