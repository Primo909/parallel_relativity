import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("Data/time.csv").drop(columns=["control"])
df=df.groupby(by=["size","N"],as_index=False).mean()
print(df)
#print(df)
#for i in range(len(df["N"])):
colors=['orange', 'blue', 'green']
for i in df["N"].unique():
    temp = df.copy().loc[df["N"]==i]
    single_core = temp.loc[df["size"]==1]["time"].values
    speedup = single_core/temp["time"]
    plt.plot(temp["size"],speedup,marker=".", label=i)
#x_lin = np.linspace(1,8,100)
#plt.plot(x_lin,x_lin,"k",label="X = 1")

plt.xlabel("# proc")
plt.ylabel("speedup")
plt.xlim(1,8)
plt.legend()
plt.show()
