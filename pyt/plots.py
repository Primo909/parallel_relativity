import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation


fig,axe= plt.subplots()

dt = 1E-4
step = 100

def animate_solution(i):
    name_file = "Data/numerical_phi_ite" + str(step*i+1)+".dat"
    iteration_data = np.loadtxt(name_file)
    x = iteration_data[:,0]
    y = iteration_data[:,1]
    axe.clear()
    axe.plot(x,y)
    axe.set_xlim([-1,1])
    axe.set_ylim([0,4])
    plt.ylabel("\u03A6")
    plt.xlabel("x")
    plt.plot(x, y, color="green")
    plt.title("T= " + str(round(dt*i*step, 2)))

ani = FuncAnimation(fig, func=animate_solution, frames=int(2/dt/step), interval=1, repeat=False)
#plt.show()
print("[ ] Building animation, estimated time: 30s")
ani.save('pyt/Figures/solution.gif', writer='Pillow', fps=100)
print("[X] Animation built")
print("Animation can be found in pyt/Figures/solution.gif")




