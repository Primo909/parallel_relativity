import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


conv_data = np.loadtxt('Data/normConvTest.dat')
t = conv_data[:,0] 
p = conv_data[:,1] 
plt.rc('axes', labelsize=16)
fig, ax = plt.subplots()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.ylim(3.985, 4.015)
plt.scatter(t, p, s=1)
plt.xlabel('t')
plt.ylabel('p')
plt.savefig("pyt/Figures/normconv_all.png")
plt.show()


fig, ax = plt.subplots()
plt.scatter(t, p, s=1)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(0.8, 1.2)
plt.ylim(3.985, 4.015)
plt.xlabel('t')
plt.ylabel('p')
plt.savefig("pyt/Figures/normconv_zoom1.png")
plt.show()
