
"""
Restricted three-body problem in the rotating frame.
Test particle trajectory in a binary star system.
"""

import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Constants
G = 6.67430e-11 # m^3 kg^-1 s^-2
M1 = 1.989e30 # mass of the first body (e.g. the Sun) in kg
M2 = 1.972e30 # mass of the second body (e.g. the Earth) in kg

# distances of two bodies from the center of mass
d = 1.496e11 * 2
M_total = M1 + M2
x1 = -M2 / M_total * d
x2 =  M1 / M_total * d
y1, y2 = 0, 0
omega = np.sqrt(G * M_total / d**3)

# initial conditions
T = 3 * (2 * np.pi / omega) # un período orbital
x0  = 1e10
vx0 = 0.0
y0  = 1e10
vy0 = 0.0


# Equations of motion in the rotating frame
def eq(t,state):
    x     = state[0]
    x_dot = state[1]
    y     = state[2]
    y_dot = state[3]
    
    x_ddot = -(G*M1*(x-x1)/((x-x1)**2 + (y-y1)**2)**(3/2))  - (G*M2*(x-x2)/((x-x2)**2 + (y-y2)**2)**(3/2)) + omega**2 * x + 2*omega*y_dot # depende de y
    y_ddot = -(G*M1*(y-y1)/((x-x1)**2 + (y-y1)**2)**(3/2))  - (G*M2*(y-y2)/((x-x2)**2 + (y-y2)**2)**(3/2)) + omega**2 * y  - 2*omega*x_dot # depende de x
    
    return [x_dot, x_ddot, y_dot, y_ddot]

sol = solve_ivp(eq, [0, T], [x0, vx0, y0, vy0], t_eval=np.linspace(0, T, 1000))

# Extracting the trajectory
x = sol.y[0]
y = sol.y[2]

# Plotting the trajectory

plt.figure(figsize=(8, 8))
plt.plot(x, y, lw=0.8, color='steelblue', label='Trajectory')  
plt.plot(x1, y1, '.', color='orange', markersize=15, label='M1')
plt.plot(x2, y2, '.', color='red', markersize=15, label='M2')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.legend()
plt.axis('equal')   
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()