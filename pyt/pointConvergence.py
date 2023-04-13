import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

title = 'L1 per cell with dt = 5e-5'
fig_title = title
conv_data = np.loadtxt('Data/gauss_zero.dat')
axis = conv_data[:,0] # the values of dx
phi_low_mid = conv_data[:,1] # the values of dx
phi_mid_high = conv_data[:,2] # the corresponding values of L1-norm

fig, ax = plt.subplots()
#majors = [1e-5, 5e-5] # used to fix the ticks in the plot
#ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
#ax.yaxis.set_major_locator(plt.MaxNLocator(3))
#plt.yscale('log')
#plt.scatter(dxs, l1s, s=1)
plt.plot(axis,phi_low_mid/16,label="phi_low_mid")
plt.plot(axis,phi_mid_high,label="phi_mid_high")
plt.title(title)
plt.xlabel('space steps')
plt.ylabel('phi differences')
plt.legend()
#plt.savefig(fig_title + '.jpg')
plt.show()
ratio = phi_low_mid/phi_mid_high
plt.plot(ratio,marker='.')
plt.ylim(0,20)
#plt.show()
plt.savefig("pyt/Figures/PointConv.pdf")
#print(np.argwhere(ratio<14))

