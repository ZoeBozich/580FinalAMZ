import matplotlib.pyplot as plt
import numpy as np
import os

from matplotlib import cm
from matplotlib.ticker import LinearLocator

from scipy.special import exp1
import sys

def generate_unique_filename(base_name, extension):
    counter = 0
    while True:
        filename = f"{base_name}_{counter}.{extension}"
        if not os.path.exists(filename):
            return filename
        counter += 1

KEY = "#################################"
DEFAULT_OUT = "output.dat"

print("Welcome to our associated python program for the final project!")

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
rho_0 = float(ins[4])
N = int(ins[5])
W = int(ins[6])

# Now that we have all of these parameters, can make properly labeled plots!
data = np.genfromtxt(DEFAULT_OUT, dtype='float', delimiter=W, skip_header=3)
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

axes.xaxis.set_label_text("x (units)")
axes.yaxis.set_label_text("y (units)")
axes.zaxis.set_label_text(r"$\rho$ (units)")
axes.set_title(fr"$\rho(x,y)$ w/ $\rho_0=${rho_0}, $R_x=${R_x}, $R_y=${R_y}, $L_x=${L_x}, $L_y=${L_y}.")

fig.colorbar(surf, shrink=0.5, aspect=5, location='left') # Maybe not needed
#Saving
rho_filename = generate_unique_filename("rho_plot", "png")
fig.savefig(rho_filename)
print(f"Charge density plot saved as {rho_filename}")

#POTENTIAL
fig, axes = plt.subplots(subplot_kw={"projection": "3d"})
surf = axes.plot_surface(x_grid, y_grid, V_grid, cmap=cm.coolwarm, linewidth=0, antialiased=False)

axes.zaxis.set_major_locator(LinearLocator(10))

axes.zaxis.set_major_formatter('{x:.02f}')

axes.xaxis.set_label_text("x (units)")
axes.yaxis.set_label_text("y (units)")
axes.zaxis.set_label_text(r"$V$ (units)")
axes.set_title(fr"$V(x,y)$ w/ $\rho_0=${rho_0}, $R_x=${R_x}, $R_y=${R_y}, $L_x=${L_x}, $L_y=${L_y}.")

fig.colorbar(surf, shrink=0.5, aspect=5, location='left') # Maybe not needed
#Saving
potential_filename = generate_unique_filename("potential_plot", "png")
fig.savefig(potential_filename)
print(f"Potential plot saved as {potential_filename}")

print(f"The maximum potential is {np.max(V)}")