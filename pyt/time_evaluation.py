# here, serial execution time is considered as single-core run of the parallelized code

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

df = pd.read_csv("Data/time.csv").drop(columns=["control"])
df =  df.groupby(by=["size","N"],as_index=False).mean()
speedup_array = []
fig, ax = plt.subplots()
print(df)


# calculating speedup, eficiency for each combination of (N, size)

plt.xlabel("# cores")
plt.ylabel("Speedup")
for i in df["N"].unique():
		temp = df.copy().loc[df["N"]==i]
		single_core = temp.loc[df["size"]==1]["time"].values
		speedup = single_core/temp["time"]
		plt.plot(temp["size"], speedup, marker=".", label="N = " + str(i))
plt.legend()
plt.show()
fig, ax = plt.subplots()
plt.xlabel("# cores")
plt.ylabel("Eficiency")

for i in df["N"].unique():
		temp = df.copy().loc[df["N"]==i]
		#size = temp["size"]
		single_core = temp.loc[df["size"]==1]["time"].values
		speedup = single_core/temp["time"]
		plt.plot(temp["size"], speedup/temp["size"], marker=".", label="N = " + str(i))

#print(np.array(speedup_array))
# line of X = 1
#x = np.linspace(1,4,100)
#plt.plot(x,x,"k",label="X=1")

# labels and limits
ax.xaxis.set_major_locator(MaxNLocator(integer=True)) # x-axis with integer values only
#plt.ylim(1,2.8)
plt.legend()
plt.savefig('pyt/Figures/speedup.pdf')
plt.show()
