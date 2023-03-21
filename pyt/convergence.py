import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

conv_data = np.loadtxt('../conv_file.dat')
dxs = conv_data[:,0] # the values of dx
l1s = conv_data[:,1] # the corresponding values of L1-norm

fig, ax = plt.subplots()
majors = [1e-5, 5e-5]
ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
ax.yaxis.set_major_locator(plt.MaxNLocator(3))
plt.scatter(dxs, l1s, s=1)
plt.xlabel('dx')
plt.ylabel('L1 - norm')
plt.show()
