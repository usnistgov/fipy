#!/usr/bin/env python

from fipy.viewers.numpyGistViewer import NumpyGistViewer
from parameters import parameters
import fipy.tools.dump as dump
from fipy.viewers.gnuplotViewer import GnuplotViewer

import Numeric

experimentalParameters = parameters['experimental parameters']

ions, acc, dis = dump.read('dump1')

ny, nx = ions.shape

symions = Numeric.zeros((ny, 2 * nx), 'd')
symacc = Numeric.zeros((ny, 2 * nx), 'd')
symdis = Numeric.zeros((ny, 2 * nx), 'd')

symions[:,:nx] = ions
symions[:,nx:] = ions[:,::-1]
symacc[:,:nx] = acc
symacc[:,nx:] = acc[:,::-1]
symdis[:,:nx] = dis
symdis[:,nx:] = dis[:,::-1]


viewers = (
    NumpyGistViewer(symions, minVal = 0., maxVal = experimentalParameters['bulk metal ion concentration'], palette = 'rainbow.gp', grid = 0, limits = (0, ny, 0, ny), dpi = 100),
    NumpyGistViewer(symacc, minVal = 0.5, maxVal = 0., palette = 'heat.gp', grid = 0, limits = (0, ny, 0, ny), dpi = 100),
    NumpyGistViewer(symdis, minVal = -1e-8, maxVal = 1e-8, palette = 'heat.gp', grid = 0, limits = (0, ny, 0, ny), dpi = 100)
    )


ionViewer = GnuplotViewer(symions, maxVal = experimentalParameters['bulk metal ion concentration'], minVal = -50. )
accViewer = GnuplotViewer(symacc, maxVal = .5, minVal = 0., palette = 'gray negative')
disViewer = GnuplotViewer(symdis, maxVal = 1e-8, minVal = -1e-8, palette = 'gray')
     
for step in range(1, 900, 20):

     nions, nacc, ndis= dump.read('dump' + str(step))

     ndis = ndis + (dis < 0) * 100.

     nions = Numeric.where(dis < 0,
                           250.,
                           Numeric.where(ndis < 0,
                                         -50.,
                                         nions))

     nacc = Numeric.where(dis < 0,
                          .5,
                          nacc)

     symions[:,:nx] = nions
     symions[:,nx:] = nions[:,::-1]
     symacc[:,:nx] = nacc
     symacc[:,nx:] = nacc[:,::-1]
     symdis[:,:nx] = ndis
     symdis[:,nx:] = ndis[:,::-1]

     ionViewer.plot(fileName = 'ion' + str(step) + '.pdf')
     accViewer.plot(fileName = 'acc' + str(step) + '.pdf')
     disViewer.plot(fileName = 'dis' + str(step) + '.pdf')

     print step

raw_input('finished')
##print ions
