#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "depositionRateVariable.py"
 #                                    created: 8/26/04 {10:39:23 AM} 
 #                                last update: 8/26/04 {4:00:40 PM} 
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
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

"""

The deposition rate variable is the distance the interface moves in a
unit of time. It is proportional to the current density such that,

.. raw:: latex

    $$ v = \\frac{i \\Omega}{n F} $$

The density is given by,

.. raw:: latex

    $$ i = i_0 \\frac{c_c^i}{c_c^{\\infty}} \\exp{ \\left( \\frac{- \\alpha F}{R T} \\eta \\right) } $$

The exchange current density is an empirical function of accelerator coverage,

.. raw:: latex

    $$ i_0(\\theta) = b_0 + b_1 \\theta $$

Here is a small test,

   >>> parameters = {
   ...     'experimental parameters' :
   ...         {
   ...             'transfer coefficient'     : 0.5,
   ...             'temperature'              : 298.,
   ...             'overpotential'            : -0.3,
   ...             'bulk metal ion concentration' : 278.,
   ...             'exchange current density' :
   ...                 {
   ...                     'constant'               : 0.26,
   ...                     'accelerator dependence' : 45.
   ...                 }
   ...         },
   ...     'material properties' :
   ...         {
   ...             'Faradays constant' : 9.6e4,
   ...             'gas constant'      : 8.314
   ...         },
   ...     'metal ion properties' :
   ...         {
   ...             'atomic volume' : 7.1e-6,
   ...             'ion charge'    : 2.
   ...         }
   ...     }
   >>> from fipy.meshes.grid2D import Grid2D
   >>> mesh = Grid2D(nx = 2, ny = 1, dx = 1., dy = 1.)
   >>> from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
   >>> distanceVariable = DistanceVariable(mesh = mesh, value = (-0.5, 0.5))
   >>> from fipy.models.levelSet.surfactant.surfactantVariable import SurfactantVariable
   >>> acceleratorVariable = SurfactantVariable(distanceVar = distanceVariable, value = (0, 1))
   >>> print DepositionRateVariable(metalIonVariable = 278., acceleratorVariable = acceleratorVariable, parameters = parameters)
   [  3.21448659e-09,  5.59567934e-07,]
   >>> print DepositionRateVariable(metalIonVariable = parameters['experimental parameters']['bulk metal ion concentration'] / 2, acceleratorVariable = acceleratorVariable, parameters = parameters)
   [  1.60724329e-09,  2.79783967e-07,]
   
"""

__docformat__ = 'restructuredtext'

import Numeric

from fipy.variables.cellVariable import CellVariable

class DepositionRateVariable(CellVariable):

    def __init__(self, metalIonVariable, acceleratorVariable = None, parameters = None, name = 'deposition rate variable'):
        """

        The following arguments are required to instatiate a
        `MetalIonSourceVariable`,

        `ionVar` - A `CellVariable`.

        `parameters` - A dictionary with the correct form as given in
        the test.

        `acceleratorCoverage` - the accelerator coverage close to the
        interface, must a `SurfactantVariable`

        """
        
        CellVariable.__init__(self, acceleratorVariable.getMesh(), hasOld = 0, name = name)

        self.metalIonVariable = self.requires(metalIonVariable)

        self.acceleratorVariable = self.requires(acceleratorVariable)
        self.parameters = parameters

        faradaysConstant = self.parameters['material properties']['Faradays constant']
        transferCoefficient = self.parameters['experimental parameters']['transfer coefficient']
        temperature = self.parameters['experimental parameters']['temperature']
        overpotential = self.parameters['experimental parameters']['overpotential']
        gasConstant = self.parameters['material properties']['gas constant']

        self.expo = Numeric.exp(- transferCoefficient * faradaysConstant * overpotential / gasConstant / temperature)

    def _calcValue(self):

        bulkConcentration = self.parameters['experimental parameters']['bulk metal ion concentration']
        atomicVolume = self.parameters['metal ion properties']['atomic volume']
        charge = self.parameters['metal ion properties']['ion charge']
        b0 = self.parameters['experimental parameters']['exchange current density']['constant']
        b1 = self.parameters['experimental parameters']['exchange current density']['accelerator dependence']
        faradaysConstant = self.parameters['material properties']['Faradays constant']

        exchangeCurrentDensity = b0 + b1 * Numeric.array(self.acceleratorVariable.getInterfaceValue())

        currentDensity = exchangeCurrentDensity * (Numeric.array(self.metalIonVariable) / bulkConcentration) * self.expo

        argmax = Numeric.argmax(currentDensity * atomicVolume / charge / faradaysConstant)
##        print "max velocity",(currentDensity * atomicVolume / charge / faradaysConstant)[argmax]
##        print "copper at max velocity",Numeric.array(self.metalIonVariable)[argmax]
##        print "length of ion variable",len(Numeric.array(self.metalIonVariable))
        self.value = currentDensity * atomicVolume / charge / faradaysConstant

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__": 
    _test()
