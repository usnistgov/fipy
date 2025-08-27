import numpy as np
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.scale import SymmetricalLogTransform
from matplotlib.ticker import (SymmetricalLogLocator, LogFormatter,
                               LogFormatterMathtext)

def plot_sparse(coo, rhs, fig, title="sparse linear system", log=False):
    fig.suptitle(title)
    sparse_grid = fig.add_gridspec(16, 24)

    L_ax = fig.add_subplot(sparse_grid[:, :16])
    c_ax = fig.add_subplot(sparse_grid[:, -5])
    b_ax = fig.add_subplot(sparse_grid[:, 16:-5])

    L = coo.todense()
    b = rhs[:, np.newaxis]

    if log:
        t = SymmetricalLogTransform(base=10, linthresh=0.0001, linscale=10)
        L = t.transform(L.flat).reshape((32, 32))
        b = t.transform(b).reshape((32, 1))
        loc = SymmetricalLogLocator(transform=t)
        fmt = LogFormatterMathtext(linthresh=0.0001)
    else:
        loc = None
        fmt = None

    zRange = max(abs(L).max(), abs(b).max())

    L_ax.matshow(L, cmap="RdBu", vmin=-zRange, vmax=zRange)
    b_ax.imshow(b, cmap="RdBu", vmin=-zRange, vmax=zRange, aspect=0.5)
    
    norm = Normalize(vmin=-zRange, vmax=zRange)
    ColorbarBase(ax=c_ax, cmap="RdBu", norm=norm, orientation='vertical',
                 format=fmt, ticks=loc)
                 
    L_ax.text(0.5, -0.1, "L",
              transform=L_ax.transAxes, horizontalalignment='center', verticalalignment='baseline')
    b_ax.text(0.5, -0.1, "b",
              transform=b_ax.transAxes, horizontalalignment='center', verticalalignment='baseline')

    for ax in [L_ax, b_ax]:
        ax.set(xticks=[], yticks=[])
