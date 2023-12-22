import matplotlib.pyplot as plt
import numpy as np

from matplotlib import cm
from matplotlib.ticker import LinearLocator

from scipy.special import exp1
import sys

KEY = "#################################"
DEFAULT_OUT = "output.dat"

print("Welcome to our associated python program for the final project potential solver!")

f = open(DEFAULT_OUT)
line1 = f.readline()
line2 = f.readline()
line3 = f.readline()
            
if not (KEY in line1 and KEY in line3):
    print("Error: incorrect data file")
    sys.exit()
    
ins = line2.split(',')
L_x = float(ins[0])
L_y = float(ins[1])
R_x = float(ins[2])
R_y = float(ins[3])
x_0 = float(ins[4])
y_0 = float(ins[5])
rho_0 = float(ins[6])
N = int(ins[7])
N_mn = int(ins[8])
M = int(ins[9])
W = int(ins[10])

# Now that we have all of these parameters, can make properly labeled plots!
data = np.genfromtxt(DEFAULT_OUT, dtype='float', delimiter=W, skip_header=4)
transposed_data = np.transpose(data)

x = transposed_data[0]
y = transposed_data[1]
rho = transposed_data[2]
V = transposed_data[3]

x_grid = np.reshape(x, (N,N))
y_grid = np.reshape(y, (N,N))
rho_grid = np.reshape(rho, (N,N))
V_grid = np.reshape(V, (N,N))

# Charge density
fig, axes = plt.subplots(subplot_kw={"projection": "3d"})
surf = axes.plot_surface(x_grid, y_grid, rho_grid, cmap=cm.coolwarm, linewidth=0, antialiased=False)

axes.zaxis.set_major_locator(LinearLocator(10))

axes.zaxis.set_major_formatter('{x:.02f}')

axes.xaxis.set_label_text(r"$x$ (m)")
axes.yaxis.set_label_text(r"$y$ (m)")
axes.zaxis.set_label_text(r"$\rho$ (C/m$^3$)")
axes.set_title(fr"$\rho(x,y)$ with $\rho_0=${rho_0}, $R_x=${R_x}, $R_y=${R_y}, $x_0=${x_0} and $y_0=${y_0}")

fig.colorbar(surf, shrink=0.5, aspect=5, location='left') # Maybe not needed
fig.savefig("rho_plot.png")

# Potential
fig, axes = plt.subplots(subplot_kw={"projection": "3d"})
surf = axes.plot_surface(x_grid, y_grid, V_grid, cmap=cm.coolwarm, linewidth=0, antialiased=False)

axes.zaxis.set_major_locator(LinearLocator(10))

axes.zaxis.set_major_formatter('{x:.02f}')

axes.xaxis.set_label_text(r"$x$ (m)")
axes.yaxis.set_label_text(r"$y$ (m)")
axes.zaxis.set_label_text(r"$V$ (C/m)")
axes.set_title(fr"$V(x,y)$ with $\rho_0=${rho_0}, $R_x=${R_x}, $R_y=${R_y}, $N={N}$, $x_0=${x_0}, $y_0=${y_0}, $N_{{mn}}={N_mn}, and $M={M}$.")

fig.colorbar(surf, shrink=0.5, aspect=5, location='left') # Maybe not needed
fig.savefig("potential_plot.png")

print(f"The maximum potential is {np.max(V)} C/m")