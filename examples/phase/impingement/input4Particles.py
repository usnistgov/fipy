#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 1/16/04 {12:00:06 PM}
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

from __future__ import nested_scopes
from fivol.examples.phase.examples.impingement.input import ImpingementSystem
import Numeric

class System4Particles(ImpingementSystem):
    def __init__(self, nx= 20, ny = 20):
        def make_circle(a,b,r,val):
            for i in range(n):
                for j in range(n):
                    ydis=(i+0.5)*dx
                    xdis=(j+0.5)*dx
                    if ((xdis-a)**2+(ydis-b)**2)<(r)**2:
                        th.u[i,j]=val
                        pf.u[i,j]=1.

        
        def circle(cell, a = 0., b = 0., r = 1.):
            x = cell.getCenter()[0]
            y = cell.getCenter()[1]
            if ((x - a)**2 + (y - b)**2) < r**2:
                return 1.
        
        def bottomRightCells(cell, Lx = 1., Ly = 1., circle = circle):
            return circle(cell, a = Lx, b = 0., r = Lx / 2.)

        def bottomLeftCells(cell, Lx = 1., Ly = 1., circle = circle):
            return circle(cell, a = 0., b = 0., r = Lx / 2.)
        
        def topRightCells(cell, Lx = 1., Ly = 1., circle = circle):
            return circle(cell, a = Lx, b = Ly, r = Lx / 2.)
        
        def topLeftCells(cell, Lx = 1., Ly = 1., circle = circle):
            return circle(cell, a = 0., b = Ly, r = Lx / 2.)

        def getAllCells(cell, Lx = 1., Ly = 1.):
            return 1.

        pi = Numeric.pi
    
        initialConditions = (
            { 'phase value' : 0., 'theta value' : -pi + 0.0001,        'func' : getAllCells },
            { 'phase value' : 1., 'theta value' : 2. * pi / 3.,        'func' : bottomLeftCells },
            { 'phase value' : 1., 'theta value' : -2. * pi / 3.,       'func' : bottomRightCells },
            { 'phase value' : 1., 'theta value' : -2. * pi / 3. + 0.3, 'func' : topLeftCells },
            { 'phase value' : 1., 'theta value' : 2. * pi / 3.,        'func' : topRightCells }
            )

        ImpingementSystem.__init__(self, nx = nx, ny = ny, initialConditions = initialConditions, steps = 10, drivingForce = 10.)

if __name__ == '__main__':
    system = System4Particles(nx = 20, ny = 20)
    system.run()
    raw_input()


