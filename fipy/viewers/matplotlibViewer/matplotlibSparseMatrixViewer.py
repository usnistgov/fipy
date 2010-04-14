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

import math

import pylab
# from pylab import ticker
from matplotlib import ticker
from scipy.io import mmio

from fipy.tools.numerix import arange, array, compress, isnan, log, log10, nan, nanmax, nanmin, sign, where, zeros, maximum, minimum

class SignedLogFormatter(ticker.LogFormatter):
    """
    Format values for log axis;

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
        sgn = sign(x)
        x = abs(x)
        x += self.threshold
        isDecade = self.is_decade(x)
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
    """
    Determine the tick locations for "log" axes that express both positive and
    negative values
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
            numdec = math.floor(limit+self.threshold)-math.ceil(self.threshold)

            if self._subs is None: # autosub
                if numdec>10: subs = array([1.0])
                elif numdec>6: subs = arange(2.0, b, 2.0)
                else: subs = arange(2.0, b)
                subs = log(subs) / log(b)
            else:
                subs = self._subs
                if numdec == 0 and len(subs) == 1:
                    subs = array(list(subs) + list(log(arange(2.0, b)) / log(b)))

            stride = 1
            while numdec/stride+1 > self.numticks:
                stride += 1

            for decadeStart in arange(math.floor(self.threshold), math.ceil(limit + self.threshold)+stride, stride):
                ticks = subs + decadeStart - self.threshold
                ticklocs.extend( sgn * ticks.compress(ticks > 0) )

        return array(ticklocs)

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

        exponent, remainder = divmod(math.log10(vmax - vmin), 1)

        if remainder < 0.5:
            exponent -= 1
        scale = 10**(-exponent)
        vmin = math.floor(scale*vmin)/scale
        vmax = math.ceil(scale*vmax)/scale

        return nonsingular(vmin, vmax)
                         
class MatplotlibSparseMatrixViewer:
    def __init__(self, title="Sparsity"):
        self.title = title
        
        self.margin = 0.1
        self.width = 0.8
        self.aspect = 1.3

        pylab.ion()
        
        fig = pylab.figure(figsize=[pylab.rcParams['figure.figsize'][0] * self.aspect, pylab.rcParams['figure.figsize'][1]])
        self.id = fig.number
        
        pylab.title(self.title)
        
    def plot(self, matrix, RHSvector, log='auto'):
        import tempfile
        import os
        
        (f, mtxName) = tempfile.mkstemp(suffix='.mtx')
        matrix.exportMmf(mtxName)
        mtx = mmio.mmread(mtxName)
##         f.close()
        os.remove(mtxName)
        
        pylab.ion()
        
        c = mtx.tocoo()
        y = c.row
        x = c.col
        z = c.data
        
        b = RHSvector
        
        if len(z) == 0:
            y = zeros((1,))
            x = zeros((1,))
            z = zeros((1,))

        zPlus = where(z > 0, log10(z), nan)
        zMinus = where(z < 0, log10(-z), nan)
        bPlus = where(b > 0, log10(b), nan)
        bMinus = where(b < 0, log10(-b), nan)

        if (log == True
            or (log == 'auto' 
                and (max(zPlus) - min(zPlus) > 2
                     or max(zMinus) - min(zMinus) > 2
                     or max(bPlus) - min(bPlus) > 2
                     or max(bMinus) - min(bMinus) > 2))):
            log = True
        else:
            log = False
            
        if log:
            zMin = nanmin((nanmin(zPlus), nanmin(zMinus), nanmin(bPlus), nanmin(bMinus)))
            zMax = nanmax((nanmax(zPlus), nanmax(zMinus), nanmax(bPlus), nanmax(bMinus)))
##             zThreshold = 0.5 # (zMax - zMin) / 5.
            
            zMin -= 0.5
            
            numdec = math.floor(zMax)-math.ceil(zMin)
            if numdec < 0:
                zMax += 0.5
            
            zPlus -= zMin
            zMinus -= zMin
            bPlus -= zMin
            bMinus -= zMin
            zRange = zMax - zMin
            
            if zRange == 0:
                zRange = nanmax(zPlus) + 1

            z = where(z > 0, zPlus, -zMinus)
            z = where(isnan(z), 0., z)
            b = where(b > 0, bPlus, -bMinus)
            b = where(isnan(b), 0., b)

            fmt = SignedLogFormatter(threshold=zMin)
            loc = SignedLogLocator(threshold=zMin)
            
        else:
            zRange = max(max(abs(z)), max(abs(b)))
        
            if zRange == 0:
                zRange = 1

            fmt = None
            loc = None
            

        N = matrix._getShape()[0]
        saveSize = pylab.rcParams['figure.figsize']
        size = pylab.rcParams['figure.dpi'] **2 * saveSize[0] * saveSize[1] / N**2

        pylab.ioff()
        
        pylab.figure(self.id)
        pylab.clf()

        pylab.delaxes()
        ax1 = pylab.axes([self.margin, self.margin, self.width, self.width])
        
        Mscat = pylab.scatter(x, y, c=z, 
                              vmin=-zRange, vmax=zRange, edgecolors='none', 
                              cmap=pylab.get_cmap('RdBu'), marker='s', s=size)
                             
        ax2 = pylab.axes([self.width + self.margin, self.margin, (self.width / self.aspect) / N, self.width], 
                         sharey=ax1)

        bscat = pylab.scatter(zeros((N,)), arange(N), c=b, 
                              vmin=-zRange, vmax=zRange, edgecolors='none', 
                              cmap=pylab.get_cmap('RdBu'), marker='s', s=size)

        pylab.setp((ax2.get_xticklabels(),
                    ax2.get_yticklabels(),
                    ax2.get_xticklines(),
                    ax2.get_yticklines()), visible=False)
        
        pylab.axes(ax1)
        pylab.axis([-0.5, N - 0.5, N - 0.5, -0.5])

        pylab.colorbar(format=fmt, ticks=loc)

        pylab.title(self.title)

        pylab.draw()
