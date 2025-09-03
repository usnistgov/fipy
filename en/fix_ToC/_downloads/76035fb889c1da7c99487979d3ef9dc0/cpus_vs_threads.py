from __future__ import unicode_literals
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np

df = pd.read_csv("threading.csv", comment='#', index_col="label")

df = df[(df['totaltime'] == 8.0) & (df['nx'] == 400)]

df = df[df['nthreads'] * df['ncpus'] == df['nslots']]

fig, ax = plt.subplots(figsize=(10,8))

colors = dict(petsc='blue', trilinos='red', pysparse='orange', scipy='green')
markers = dict(petsc='x', trilinos='+', pysparse='^', scipy='v')
linestyles = {1: '-', 2: '--', 4: '-.', 16: ':'}

for solver, group1 in df.groupby('solver'):
    stats = group1.groupby('nthreads')
    
    group1mean = stats.get_group(1).solvetime.mean(numeric_only=True)
    av = group1mean / stats.solvetime.mean(numeric_only=True)
    st = group1mean / stats.solvetime.std()
    st[st.isna()] = 0.
    
    ax.errorbar(av.index, av, # yerr=st, 
                marker=markers[solver], color=colors[solver], linewidth=2,
                markersize=12, label=solver)
                    
plt.xscale('log')
plt.yscale('log')

for axis in [ax.xaxis, ax.yaxis]:
    axis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))

plt.legend(loc="upper right", frameon=False, handlelength=5)

plt.xlabel("# threads")
plt.ylabel(r'"speedup" ($t_1 / t_N$)')

plt.show()
