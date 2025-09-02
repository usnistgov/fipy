import numpy as np
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.scale import SymmetricalLogTransform
from matplotlib.ticker import (SymmetricalLogLocator, LogFormatter,
                               LogFormatterMathtext)

class SparsePlot:
    def __init__(self, coo, rhs, fig, title="sparse linear system"):
        self.coo = coo
        self.rhs = rhs
        self.fig = fig
        self.title = title

        self.sparse_grid = self.fig.add_gridspec(16, 24)
        self.L_ax = self.fig.add_subplot(self.sparse_grid[:, :16])
        self.c_ax = self.fig.add_subplot(self.sparse_grid[:, -5])
        self.b_ax = self.fig.add_subplot(self.sparse_grid[:, 16:-5])

    def transform_Lb(self, L, b):
        return L, b, None, None

    def plot(self):
        self.fig.suptitle(self.title)

        L = self.coo.todense()
        b = self.rhs[:, np.newaxis]

        L, b, loc, fmt = self.transform_Lb(L, b)

        zRange = max(abs(L).max(), abs(b).max())

        self.L_ax.matshow(L, cmap="RdBu", vmin=-zRange, vmax=zRange)
        self.b_ax.imshow(b, cmap="RdBu", vmin=-zRange, vmax=zRange, aspect=0.5)

        norm = Normalize(vmin=-zRange, vmax=zRange)
        ColorbarBase(ax=self.c_ax, cmap="RdBu", norm=norm, orientation='vertical',
                     format=fmt, ticks=loc)

        self.L_ax.text(0.5, -0.1, "L",
                       transform=self.L_ax.transAxes,
                       horizontalalignment='center',
                       verticalalignment='baseline')
        self.b_ax.text(0.5, -0.1, "b",
                       transform=self.b_ax.transAxes,
                       horizontalalignment='center',
                       verticalalignment='baseline')

        for ax in [self.L_ax, self.b_ax]:
            ax.set(xticks=[], yticks=[])

class LogSparsePlot(SparsePlot):
    def transform_Lb(self, L, b):
        t = SymmetricalLogTransform(base=10, linthresh=0.0001, linscale=10)
        L = t.transform(L.flat).reshape((32, 32))
        b = t.transform(b).reshape((32, 1))
        loc = SymmetricalLogLocator(transform=t)
        fmt = LogFormatterMathtext(linthresh=0.0001)

        return L, b, loc, fmt
