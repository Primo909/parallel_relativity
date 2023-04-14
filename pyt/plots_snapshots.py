import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import sqrt

dt = 1E-4
step = 100
sigma = 0.1


list_times = [1, 50, 100, 150]
plt.rc('axes', labelsize=16)
fig, ax = plt.subplots(2, 2)
j = -1

for i in list_times:
    j += 1
    data = np.loadtxt("Data/numerical_phi_ite" + str(step*i + 1) + ".dat")
    subplot_number = "ax" + str(j)
    time = data[:,0] 
    phi = data[:,1] 

    if j == 0: 
        x=0
        y=0
    if j == 1: 
        x=1
        y=0
    if j == 2: 
        x=0
        y=1
    if j == 3: 
        x=1
        y=1

    ax[x, y].plot(time, phi, color='black')
    ax[x, y].set_ylim(0, 4)
    ax[x, y].set_title("t = " + str(0.01*i))
    ax[x, y].set_xticks([])
    ax[x, y].set_yticks([])

    if j == 1 or j == 3: 
        ax[x, y].set_xticks([-1.,-0.5, 0, 0.5, 1.])
        ax[x, y].set_xlabel('x')

    if j == 0 or j == 1:
        ax[x, y].set_yticks([0.,1., 2., 3., 4.])
        ax[x, y].set_ylabel('$\Phi$')

plt.tight_layout()
plt.savefig("pyt/Figures/snapshots.png")
plt.show()

analyt = 