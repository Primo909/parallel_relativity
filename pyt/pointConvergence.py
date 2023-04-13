import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


conv_data = np.loadtxt('Data/pointConvTest.dat')
axis = conv_data[:,0] 
phi_low_mid = conv_data[:,1] 
phi_mid_high = conv_data[:,2] 

plt.rc('axes', labelsize=16)
fig, ax = plt.subplots()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.plot(axis, phi_low_mid/16, label="$\Phi_{mid} - \Phi_{low}$")
plt.plot(axis, phi_mid_high, label="$\Phi_{mid} - \Phi_{high}$")
plt.xlabel('$x$')
plt.ylabel('$\Phi$ differences')
plt.legend()
plt.savefig('pyt/Figures/point_conv.png')
plt.show()

