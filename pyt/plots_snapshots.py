import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


dt = 1E-4
step = 100

list_times = [1, 50, 100, 150]
plt.rc('axes', labelsize=16)


for i in list_times:

    data = np.loadtxt("Data/numerical_phi_ite" + str(step*i + 1) + ".dat")
    time = data[:,0] 
    phi = data[:,1] 

    fig, ax = plt.subplots()
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_xticks([-1.,-0.5, 0, 0.5, 1.])
    plt.ylim(0, 4)
    plt.xlabel('t')
    plt.ylabel('$\Phi$')
    plt.plot(time, phi, label="t = " + str(0.01*i), color='black')
    plt.legend(loc='upper right')
    plt.savefig("pyt/Figures/snapshots" + str(i) + ".png")

plt.show()