from __future__ import unicode_literals
from matplotlib import pyplot as plt
import pandas as pd

df = pd.read_csv("threadanalyze.csv", comment='#')

df = df[df['nthreads'] * df['ncpus'] == df['nslots']]

slots16 = df[df['nslots'] == 16]
slots32 = df[df['nslots'] == 32]

av = slots16.groupby('nthreads').mean()
sd = slots16.groupby('nthreads').std()

# plt.plot(slots16['nthreads'], slots16['solvetime'], linestyle="", marker="+", color='blue', label="16 slots")
plt.errorbar(av['solvetime'].index, av['solvetime'].values, yerr=sd['solvetime'].values, linestyle="", marker="o", color='blue', label="16 slots")
plt.plot(slots32['nthreads'], slots32['solvetime'], linestyle="", marker="+", color='red', label="32 slots")

plt.legend(loc="upper right")

plt.xlabel("# threads")
plt.ylabel("time / s")

plt.show()
