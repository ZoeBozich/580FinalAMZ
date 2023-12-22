import matplotlib.pyplot as plt
import numpy as np

from matplotlib import cm
from matplotlib.ticker import LinearLocator

from scipy.special import exp1
import sys

def s_func(num):
    if num == 1:
        return ""
    else:
        return "s"

KEY = "#################################"
DEFAULT_OUT = "output.dat"

print("Welcome to our associated python program for the final project thread plotter!")

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
num_Ns = int(ins[7])
N_mn = int(ins[8])
M = int(ins[9])
W = int(ins[10])

# Now that we have all of these parameters, can make properly labeled plots!
data = np.genfromtxt(DEFAULT_OUT, dtype='float', delimiter=W, skip_header=4)
transposed_data = np.transpose(data)

threads = transposed_data[0]
N = transposed_data[1]
t = transposed_data[2]

num_threads = len(threads) // num_Ns

thread_grid = np.reshape(threads, (num_threads,num_Ns))
N_grid = np.reshape(N, (num_threads,num_Ns))
t_grid = np.reshape(t, (num_threads,num_Ns))

for i in range(num_threads):
    th = int(thread_grid[i][0])
    plt.plot(N_grid[i],t_grid[i],label=f"{th} thread{s_func(th)}")

plt.xlabel(r"Grid points $N$")
plt.ylabel(r"Runtime $t$ (s)")
plt.suptitle("Dependence of computation time on grid size and thread count")
plt.title(rf"($\rho_0={rho_0}$, $R_x={R_x}$, $R_y={R_y}$, $x_0={x_0}$, $y_0={y_0}$, $N_{{mn}}={N_mn}$, and $M={M}$)")
plt.legend()
plt.savefig("N_versus_threads.png")
print("Plot saved!")