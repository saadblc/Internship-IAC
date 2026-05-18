import numpy as np
import matplotlib.pyplot as plt

G = 1

# Masas
M1 = 1.0
M2 = 0.25

# Posiciones
x1, y1 = -1, 0
x2, y2 =  1, 0

# Velocidad angular
Omega = 1

# Grid
x = np.linspace(-3, 3, 1000)
y = np.linspace(-3, 3, 1000)

X, Y = np.meshgrid(x, y)

# Regularización
eps = 1e-3

# Distancias
r1 = np.sqrt((X - x1)**2 + (Y - y1)**2 + eps)
r2 = np.sqrt((X - x2)**2 + (Y - y2)**2 + eps)

# Potencial efectivo
Phi = (
    -G * M1 / r1
    -G * M2 / r2
    -0.5 * Omega**2 * (X**2 + Y**2)
)

# Figura
plt.figure(figsize=(8,8))

# Niveles cuidadosamente elegidos
levels = np.linspace(-3, -1.2, 40)

contours = plt.contour(
    X,
    Y,
    Phi,
    levels=levels,
    cmap='plasma'
)

plt.clabel(contours, fontsize=8)

# Dibujar masas
plt.scatter([x1, x2], [y1, y2], color='red', s=100)

plt.xlabel('x')
plt.ylabel('y')

plt.axis('equal')

plt.title('Potencial efectivo y lóbulos de Roche')

plt.show()