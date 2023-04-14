# here, serial execution time is considered as single-core run of the parallelized code

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


df = pd.read_csv("Data/time.csv").drop(columns=["control"])
df =  df.groupby(by=["size","N"],as_index=False).mean()
print(df)
plt.rc('axes', labelsize=16) 
fig, ax = plt.subplots() 
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("# cores")
plt.ylabel("Speedup")

# calculating speedup, eficiency for each combination of (N, size)

plt.xlabel("# cores")
plt.ylabel("Speedup")
for i in df["N"].unique():
		temp = df.copy().loc[df["N"]==i]
		single_core = temp.loc[df["size"]==1]["time"].values
		speedup = single_core/temp["time"]
		plt.plot(temp["size"], speedup, marker=".", label="N = " + str(i))


x = np.linspace(1,4,100)
plt.plot(x, x, "k", label="$S_n$ = n")  # line of X = X
ax.xaxis.set_major_locator(MaxNLocator(integer=True)) # x-axis with integer values only
plt.ylim(1, 2.8)
plt.legend()
plt.savefig('pyt/Figures/speedup.png')
plt.show()



fig, ax = plt.subplots()
#plt.rc('axes', labelsize=16) 
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("#cores")
plt.ylabel("Efficiency")

for i in df["N"].unique():
		temp = df.copy().loc[df["N"]==i]
		#size = temp["size"]
		single_core = temp.loc[df["size"]==1]["time"].values
		speedup = single_core/temp["time"]
		plt.plot(temp["size"], speedup/temp["size"], marker=".", label="N = " + str(i))

# labels and limits
ax.xaxis.set_major_locator(MaxNLocator(integer=True)) # x-axis with integer values only
plt.legend()
plt.savefig('pyt/Figures/efficiency.png')
plt.show()
