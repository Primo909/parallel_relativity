# here, serial execution time is considered as the time to run of the unparallelized code (main.exe, with the classes etc)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("Data/time.csv").drop(columns=["control"])
df_serial = pd.read_csv("Data/time_serial.csv").drop(columns=["control"])
df=df.groupby(by=["size","N"],as_index=False).mean()
df_serial=df_serial.groupby(by=["N"],as_index=False).mean()
#print(df)
#print(df_serial)

for i in df["N"].unique():
    if i < 5000:
        temp = df.copy().loc[df["N"]==i]
        temp_serial = df_serial.copy().loc[df_serial["N"]==i]
        single_core = temp_serial["time"].values
        speedup = single_core/temp["time"]
        plt.plot(temp["size"],speedup,marker=".", label=i)


#x_lin = np.linspace(1,8,100)
#plt.plot(x_lin,x_lin,"k",label="X = 1")
plt.xlabel("# proc")
plt.ylabel("speedup")
plt.xlim(1,8)
plt.legend()
plt.savefig('pyt/speedup_serial.pdf')
plt.show()
