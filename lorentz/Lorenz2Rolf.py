

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import scipy


fig = plt.figure()
ax = Axes3D(fig)

def lorenz_int(initial, t):
  x = initial[0]
  y = initial[1]
  z = initial[2]
  #
  sigma = 10
  rho = 28
  beta = 8.0/3
  x_dot = sigma * (y - x)
  y_dot = x * (rho -z) - y
  z_dot = x * y - beta* z
  return (x_dot, y_dot, z_dot)
    
initial = [0, 2, 1.05]

t = scipy.arange(0, 100, 0.01)
lorenz_sol = odeint(lorenz_int, initial, t)
lorenz = lorenz_int(initial, t)

x = [i[0] for i in lorenz_sol]
y = [i[1] for i in lorenz_sol]
z = [i[2] for i in lorenz_sol]

# ax.xlabel('x')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

ax.plot(x, y, z, label="lorenz" )
# ax.plot(x, y, z, label="lorenz", linestyle='none', marker='.')
ax.legend()


plt.show()
