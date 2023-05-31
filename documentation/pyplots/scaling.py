from __future__ import unicode_literals
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np

df = pd.read_csv("scaling.csv", comment='#', index_col='label')

df = df[(df['totaltime'] == 8.0) & (df['nx'] == 1000)]
df = df[df['nthreads'] * df['ncpus'] == df['nslots']]

baseline = df[df['solver'] == "pysparse"].iloc[0].solvetime
df['speedup'] = baseline / df['solvetime']

fig, ax = plt.subplots(figsize=(10,8))

colors = dict(petsc='blue', trilinos='red', pysparse='orange', scipy='green')
markers = dict(petsc='x', trilinos='+', pysparse='^', scipy='v')
linestyles = {1: '-', 2: '--', 4: '-.', 16: ':'}

for solver, group1 in df.groupby('solver'):
    for nthreads, group2 in group1.groupby('nthreads'):
        stats = group2.groupby('nslots')
        speedup = stats.mean(numeric_only=True).speedup
        yerr = stats.std(numeric_only=True).speedup
        yerr[yerr.isna()] = 0.

        ax.errorbar(speedup.index, speedup, yerr=yerr,
                    marker=markers[solver], color=colors[solver], linestyle=linestyles[nthreads], linewidth=2,
                    markersize=12, label="{} - {:.0f} thread(s)".format(solver, nthreads))
                    
plt.xscale('log')
plt.yscale('log')

for axis in [ax.xaxis, ax.yaxis]:
    axis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))

plt.legend(loc="lower right", frameon=False, handlelength=5)

plt.xlabel("# tasks")
plt.ylabel(r"speedup ($t_{\mathrm{PySparse}} / t_N$)")

plt.show()
