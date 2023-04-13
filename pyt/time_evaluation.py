# here, serial execution time is considered as single-core run of the parallelized code

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("Data/time.csv").drop(columns=["control"])
df=df.groupby(by=["size","N"],as_index=False).mean()
print(df)

for i in df["N"].unique():
	    temp = df.copy().loc[df["N"]==i]
	    single_core = temp.loc[df["size"]==1]["time"].values
	    speedup = single_core/temp["time"]
	    plt.plot(temp["size"],speedup,marker=".", label=i)

x=np.linspace(1,4,100)
plt.plot(x,x,"k",label="X=1")
plt.xlabel("# proc")
plt.ylabel("speedup")
plt.ylim(1,2.5)
plt.legend()
plt.savefig('pyt/Figures/speedup.pdf')
plt.show()
