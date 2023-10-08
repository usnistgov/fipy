from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from builtins import object
from builtins import zip
from matplotlib import cm
from matplotlib import pyplot
from matplotlib import rcParams
from matplotlib import ticker
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from scipy.io import mmio

from fipy.tools import numerix

__all__ = ["MatplotlibSparseMatrixViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class MatplotlibSparseMatrixViewer(object):
    """Displays :class:`~fipy.matrices.sparseMatrix._SparseMatrix` objects using Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """
    def __init__(self, title="Sparsity"):
        self.title = title

        self.L_width = 0.8
        self.margin = (1. - self.L_width) / 2
        self.b_width = self.margin
        self.c_width = self.margin / 3
        self.buffer = 1.5 * self.margin
        self.aspect = (self.margin + self.L_width                   # M
                       + self.buffer + self.c_width                 # colorbar
                       + self.buffer + self.b_width + self.margin)  # b

        pyplot.ion()

        fig = pyplot.figure(figsize=pyplot.figaspect(1. / self.aspect))
        self.id = fig.number

    def plot(self, matrix, RHSvector, log='auto'):
        import tempfile
        import os

        if "print" in os.environ.get('FIPY_DISPLAY_MATRIX', '').lower().split():
            print("-"*75)
            print(self.title)
            print("-"*75)
            print("L:")
            print(matrix)
            print("b:", RHSvector)

        (f, mtxName) = tempfile.mkstemp(suffix='.mtx')
        matrix.exportMmf(mtxName)
        mtx = mmio.mmread(mtxName)
        os.remove(mtxName)

        pyplot.ion()

        c = mtx.tocoo()
        y = c.row
        x = c.col
        z = c.data
        N = matrix._shape[0]

        b = RHSvector
        if numerix.shape(b) == ():
            b = numerix.zeros((N,), 'l')

        if len(z) == 0:
            y = numerix.zeros((1,), 'l')
            x = numerix.zeros((1,), 'l')
            z = numerix.zeros((1,), 'l')

        def signed_to_logs(v):
            return (numerix.where(v > 0, numerix.log10(v), numerix.nan),
                    numerix.where(v < 0, numerix.log10(-v), numerix.nan))

        def logs_to_signed(v, plus, minus):
            v = numerix.where(v > 0, plus, -minus)
            v = numerix.where(numerix.isnan(v), 0., v)

            return v

        zPlus, zMinus = signed_to_logs(z)
        bPlus, bMinus = signed_to_logs(b)

        logs = (zPlus, zMinus, bPlus, bMinus)

        log = ((log == True)
               or (log == 'auto'
                   and (numerix.nanmax(numerix.concatenate(logs))
                        - numerix.nanmin(numerix.concatenate(logs)) > 2)))

        if log:
            zMin = numerix.nanmin(numerix.concatenate(logs))
            zMax = numerix.nanmax(numerix.concatenate(logs))

            zMin -= 0.5

            numdec = numerix.floor(zMax) - numerix.ceil(zMin)
            if numdec < 0:
                zMax += 0.5

            for v in logs:
                v -= zMin

            zRange = zMax - zMin

            if zRange == 0:
                zRange = numerix.nanmax(zPlus) + 1

            z = logs_to_signed(z, zPlus, zMinus)
            b = logs_to_signed(b, bPlus, bMinus)

            fmt = None
            loc = ticker.SymmetricalLogLocator(linthresh=1, base=10)

        else:
            zRange = max(abs(numerix.concatenate((z, b))))

            if zRange == 0:
                zRange = 1

            fmt = None
            loc = None


        pyplot.ioff()

        fig = pyplot.figure(self.id)
        fig.clf()

        usetex = rcParams['text.usetex']
        rcParams['text.usetex'] = False

        cmap = cm.RdBu

        norm = Normalize(vmin=-zRange, vmax=zRange)

        x0 = self.margin
        L_ax = fig.add_axes([x0 / self.aspect, self.margin, self.L_width / self.aspect, self.L_width])
        L_ax.text(0.5, -0.1, "L",
                  transform=L_ax.transAxes, horizontalalignment='center', verticalalignment='baseline')

        x0 += self.L_width + self.buffer
        c_ax = fig.add_axes([x0 / self.aspect, self.margin, self.c_width / self.aspect, self.L_width])

        x0 += self.c_width + self.buffer
        b_ax = fig.add_axes([x0 / self.aspect, self.margin, self.b_width / self.aspect, self.L_width],
                            sharey=L_ax)
        b_ax.text(0.5, -0.1, "b",
                  transform=b_ax.transAxes, horizontalalignment='center', verticalalignment='baseline')



        def scatterRectangles(x, y, z, norm=None, cmap=None):
            patches = [Rectangle(numerix.array([X - 0.5, Y - 0.5]), 1., 1.,
                                 edgecolor='none') for X, Y in zip(x, y)]

            collection = PatchCollection(patches, norm=norm, cmap=cmap,
                                         edgecolors='none')
            collection.set_array(z)

            return collection

        L_ax.add_collection(scatterRectangles(x=x, y=y, z=z,
                                              norm=norm, cmap=cmap))

        b_ax.add_collection(scatterRectangles(x=numerix.zeros((N,), 'l'), y=numerix.arange(N), z=b,
                                              norm=norm, cmap=cmap))

        ColorbarBase(ax=c_ax, cmap=cmap, norm=norm, orientation='vertical',
                     format=fmt, ticks=loc)

        pyplot.setp((b_ax.get_xticklabels(),
                     b_ax.get_yticklabels(),
                     b_ax.get_xticklines(),
                     b_ax.get_yticklines()), visible=False)

        L_ax.set_xlim(xmin=-0.5, xmax=N-0.5)
        L_ax.set_ylim(ymax=-0.5, ymin=N-0.5)

        b_ax.set_xlim(xmin=-0.5, xmax=0.5)
        b_ax.set_ylim(ymax=-0.5, ymin=N-0.5)

        fig.suptitle(self.title, x=0.5, y=0.95, fontsize=14)

        pyplot.draw()

        rcParams['text.usetex'] = usetex
