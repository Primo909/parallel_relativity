import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import sqrt

dt = 1E-4
step = 100

list_times = [1, 50, 100, 150]
#plt.rc('axes', labelsize=16)
fig = plt.figure()
ax = fig.add_subplot(111)  
ax0 = fig.add_subplot(221)
ax1 = fig.add_subplot(222)
ax2 = fig.add_subplot(223)
ax3 = fig.add_subplot(224)
j = -1

ax.tick_params(left=False)

for i in list_times:
    j += 1
    data = np.loadtxt("Data/numerical_phi_ite" + str(step*i + 1) + ".dat")
    subplot_number = "ax" + str(j)
    time = data[:,0] 
    phi = data[:,1] 

    #plt.xticks(fontsize=14)
    #plt.yticks(fontsize=14)
    plt.ylim(0, 4)
    #plt.xlabel('t')
    #plt.ylabel('$\Phi$')
    (eval(subplot_number)).plot(time, phi, label="t = " + str(0.01*i), color='black')
    (eval(subplot_number)).set_xticks([-1.,-0.5, 0, 0.5, 1.])
    plt.legend(loc='upper right')
    plt.savefig("pyt/Figures/snapshots" + str(i) + ".png")

plt.show()