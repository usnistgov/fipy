## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 #
 # FILE: "matplotlibSparseMatrixViewer.py"
 #
 # Author: Jonathan Guyer <guyer@nist.gov>
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #
 # #############################################################################
 ##

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

class SignedLogFormatter(ticker.LogFormatter):
    """Format signed values for log axis

    formatted values range from large positive to small positive to small negative to large negative.

    if attribute decadeOnly is True, only the decades will be labelled.
    """
    def __init__(self, base=10.0, labelOnlyBase=True, threshold=0.):
        """
        base is used to locate the decade tick,
        which will be the only one to be labeled if labelOnlyBase
        is False
        """
        self._base = base+0.0
        self.labelOnlyBase=labelOnlyBase
        self.threshold = threshold

    def __call__(self, x, pos=None):
        'Return the format for tick val x at position pos'
        vmin, vmax = self.axis.get_view_interval()
        d = abs(vmax - vmin)
        b=self._base
        # only label the decades
        sgn = numerix.sign(x)
        x = abs(x)
        x += self.threshold
        isDecade = ticker.is_decade(x)
        bx = b**x
        sgnbx = sgn * bx
        if not isDecade and self.labelOnlyBase: s = ''
        elif x > 4: s= '%+1.0e' % sgnbx
        elif x < 0: s =  '%+1.0e' % sgnbx
        else        : s =  self.pprint_val(sgnbx, d)
        return s

    def pprint_val(self, x, d):
        #if the number is not too big and it's an int, format it as an
        #int
        if abs(x)<1e4 and x==int(x): return '%+d' % x

        d = 10**(d / 2.)

        if d < 1e-2: fmt = '%+1.3e'
        elif d < 1e-1: fmt = '%+1.3f'
        elif d > 1e5: fmt = '%+1.1e'
        elif d > 10 : fmt = '%+1.1f'
        elif d > 1 : fmt = '%+1.2f'
        else: fmt = '%+1.3f'
        s =  fmt % x
        #print d, x, fmt, s
        tup = s.split('e')
        if len(tup)==2:
            mantissa = tup[0].rstrip('0').rstrip('.')
            sgn = tup[1][0].replace('+', '')
            exponent = tup[1][1:].lstrip('0')
            s = '%se%s%s' %(mantissa, sgn, exponent)
        else:
            s = s.rstrip('0').rstrip('.')

        return s

class SignedLogLocator(ticker.LogLocator):
    """Determine the tick locations for "signed" log axes

    Locate ticks from large positive to small positive to small negative to large negative
    """

    def __init__(self, base=10.0, subs=[1.0], threshold=0.):
        """
        place ticks on the location= base**i*subs[j]
        """
        ticker.LogLocator.__init__(self, base, subs)
        self.numticks = 7
        self.threshold = threshold

    def _set_numticks(self):
        self.numticks = 7  # todo; be smart here; this is just for dev

    def __call__(self):
        'Return the locations of the ticks'
        b=self._base

        vmin, vmax = self.axis.get_view_interval()

        if vmax<vmin:
            vmin, vmax = vmax, vmin

        if vmax * vmin > 0:
            raise ValueError("The interval must range from negative to positive values")

        ticklocs = []

        for limit, sgn in zip([vmax, -vmin], [1, -1]):
            numdec = numerix.floor(limit+self.threshold)-numerix.ceil(self.threshold)

            if self._subs is None: # autosub
                if numdec>10:
                    subs = numerix.array([1.0])
                elif numdec>6:
                    subs = numerix.arange(2.0, b, 2.0)
                else:
                    subs = numerix.arange(2.0, b)
                subs = numerix.log(subs) / numerix.log(b)
            else:
                subs = self._subs
                if numdec == 0 and len(subs) == 1:
                    subs = numerix.array(list(subs) + list(numerix.log(numerix.arange(2.0, b)) / numerix.log(b)))

            stride = 1
            while numdec/stride+1 > self.numticks:
                stride += 1

            for decadeStart in numerix.arange(numerix.floor(self.threshold),
                                              numerix.ceil(limit + self.threshold)+stride,
                                              stride):
                ticks = subs + decadeStart - self.threshold
                ticklocs.extend( sgn * ticks.compress(ticks > 0) )

        return numerix.array(ticklocs)

    def autoscale(self):
        'Try to choose the view limits intelligently'
        self.verify_intervals()
        vmin, vmax = self.dataInterval.get_bounds()

        if vmax<vmin:
            vmin, vmax = vmax, vmin

        if vmax * vmin > 0:
            raise ValueError("The interval must range from negative to positive values")

        if vmax == self.threshold:
            vmax += 1
        if vmin == -self.threshold:
            vmin -= 1

        exponent, remainder = divmod(numerix.log10(vmax - vmin), 1)

        if remainder < 0.5:
            exponent -= 1
        scale = 10**(-exponent)
        vmin = numerix.floor(scale*vmin)/scale
        vmax = numerix.ceil(scale*vmax)/scale

        return nonsingular(vmin, vmax)

class MatplotlibSparseMatrixViewer:
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

        if "print" in os.environ['FIPY_DISPLAY_MATRIX'].lower().split():
            print "-"*75
            print self.title
            print "-"*75
            print "L:"
            print matrix
            print "b:", RHSvector

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

            fmt = SignedLogFormatter(threshold=zMin)
            loc = SignedLogLocator(threshold=zMin)

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
