import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Parámetros
G = 1
m1 = 1
m2 = 1
M = m1 + m2

# Tiempo
t = np.linspace(0, 40, 1000)

# Radio decreciente (pérdida de momento angular)
a0 = 5
a = a0 * np.exp(-0.03 * t)

# Frecuencia angular creciente
omega = np.sqrt(M / a**3)

# Fase orbital
theta = np.cumsum(omega) * (t[1] - t[0])

# Posiciones
x1 = (m2/M) * a * np.cos(theta)
y1 = (m2/M) * a * np.sin(theta)

x2 = -(m1/M) * a * np.cos(theta)
y2 = -(m1/M) * a * np.sin(theta)

# Figura
fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlim(-6,6)
ax.set_ylim(-6,6)

star1, = ax.plot([], [], 'o', markersize=12)
star2, = ax.plot([], [], 'o', markersize=12)

trail1, = ax.plot([], [], lw=1)
trail2, = ax.plot([], [], lw=1)

def update(frame):
    star1.set_data([x1[frame]], [y1[frame]])
    star2.set_data([x2[frame]], [y2[frame]])

    trail1.set_data(x1[:frame], y1[:frame])
    trail2.set_data(x2[:frame], y2[:frame])

    return star1, star2, trail1, trail2

ani = FuncAnimation(fig, update, frames=len(t), interval=20)

plt.show()
