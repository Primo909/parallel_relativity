# here, serial execution time is considered as single-core run of the parallelized code

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("Data/time.csv").drop(columns=["control"])
df=df.groupby(by=["size","N"],as_index=False).mean()
print(df)

for i in df["N"].unique():
    if i < 5000:
	    temp = df.copy().loc[df["N"]==i]
	    single_core = temp.loc[df["size"]==1]["time"].values
	    speedup = single_core/temp["time"]
	    plt.plot(temp["size"],speedup,marker=".", label=i)

plt.xlabel("# proc")
plt.ylabel("speedup")
plt.xlim(1,8)
plt.legend()
plt.savefig('pyt/Figures/speedup.pdf')
plt.show()
