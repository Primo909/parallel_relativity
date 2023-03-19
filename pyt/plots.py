import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

fig,axe= plt.subplots()

def animate_solution(i):
    name_file= "../Data/numerical_phi_ite"+str(i+1)+".dat"
    iteration_data= np.loadtxt(name_file)
    x= iteration_data[:,0]
    y= iteration_data[:,1]
    axe.clear()
    axe.plot(x,y)
    axe.set_xlim([0,1])
    axe.set_ylim([0,2])

ani = FuncAnimation(fig, func=animate_solution,frames=1000,interval=1, repeat=False)

plt.show()


