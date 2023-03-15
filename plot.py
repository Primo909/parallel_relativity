# Standardpakete
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('out.csv')
df.plot()
#plt.ylim(-1,1)
plt.show()
