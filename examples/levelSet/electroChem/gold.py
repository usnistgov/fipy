#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "gold.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
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
 # ###################################################################
 ##

r"""Model electrochemical superfill of gold using the CEAC mechanism.

This input file
is a demonstration of the use of :term:`FiPy` for
modeling gold superfill. The material properties and experimental
parameters used are roughly those that have been previously
published :cite:`NIST:damascene:2005`.

To run this example from the base fipy directory type::

    $ python examples/levelSet/electroChem/gold.py

at the command line. The results of the simulation will be displayed and the word
``finished`` in the terminal at the end of the simulation. The simulation will
only run for 10 time steps. To run with a different number of time steps change
the ``numberOfSteps`` argument as follows,

.. index:: runGold

>>> runGold(numberOfSteps=10, displayViewers=False) # doctest: +GMSH, +LSM
1

Change the ``displayViewers`` argument to ``True`` if you wish to see the
results displayed on the screen. This example has a more realistic
default boundary layer depth and thus requires `gmsh` to construct a
more complex mesh.

.. index:: gmsh

There are a few differences between the gold superfill model presented
in this example and in :mod:`examples.levelSet.electroChem.simpleTrenchSystem`.
Most default
values have changed to account for a different metal ion (gold)
and catalyst (lead). In this system the catalyst is not present in
the electrolyte but instead has a non-zero initial coverage. Thus
quantities associated with bulk catalyst and catalyst accumulation
are not defined. The current density is given by,

.. math::

   i = \frac{c_m}{c_m^{\infty}} \left( b_0 + b_1 \theta \right).

The
more common representation of the current density includes an
exponential part. Here it is buried in :math:`b_0` and :math:`b_1`. The
governing equation for catalyst evolution includes a term for
catalyst consumption on the interface and is given by

.. math::

   \dot{\theta} = J v \theta - k_c v \theta

where :math:`k_c` is the consumption coefficient
(``consumptionRateConstant``). The trench geometry is also given a
slight taper, given by ``taperAngle``.

If the MayaVi plotting software is
installed (see :ref:`INSTALLATION`) then a plot should
appear that is updated every 10 time steps and will eventually
resemble the image below.

.. image:: electroChem/inputGold.*
   :width: 90%
   :align: center
   :alt: catalyst coverage as a function of time during gold superfill

.. .. bibmissing:: /documentation/refs.bib
    :sort:
"""
__docformat__ = 'restructuredtext'

from fipy import CellVariable, SurfactantVariable, TransientTerm, FirstOrderAdvectionTerm, MultiViewer, Viewer
from fipy.tools import numerix
from trenchMesh import TrenchMesh
from gapFillDistanceVariable  import GapFillDistanceVariable
from metalIonDiffusionEquation import buildMetalIonDiffusionEquation
from adsorbingSurfactantEquation import AdsorbingSurfactantEquation

def runGold(faradaysConstant=9.6e4,
            consumptionRateConstant=2.6e+6,
            molarVolume=10.21e-6,
            charge=1.0,
            metalDiffusion=1.7e-9,
            metalConcentration=20.0,
            catalystCoverage=0.15,
            currentDensity0=3e-2 * 16,
            currentDensity1=6.5e-1 * 16,
            cellSize=0.1e-7,
            trenchDepth=0.2e-6,
            aspectRatio=1.47,
            trenchSpacing=0.5e-6,
            boundaryLayerDepth=90.0e-6,
            numberOfSteps=10,
            taperAngle=6.0,
            displayViewers=True):

    cflNumber = 0.2
    numberOfCellsInNarrowBand = 20

    mesh = TrenchMesh(cellSize = cellSize,
                      trenchSpacing = trenchSpacing,
                      trenchDepth = trenchDepth,
                      boundaryLayerDepth = boundaryLayerDepth,
                      aspectRatio = aspectRatio,
                      angle = numerix.pi * taperAngle / 180.)

    narrowBandWidth = numberOfCellsInNarrowBand * cellSize

    distanceVar = GapFillDistanceVariable(
       name = 'distance variable',
       mesh = mesh,
       value = -1.)

    distanceVar.setValue(1., where=mesh.electrolyteMask)
    distanceVar.calcDistanceFunction()

    catalystVar = SurfactantVariable(
        name = "catalyst variable",
        value = catalystCoverage,
        distanceVar = distanceVar)

    metalVar = CellVariable(
        name = 'metal variable',
        mesh = mesh,
        value = metalConcentration)

    exchangeCurrentDensity = currentDensity0 + currentDensity1 * catalystVar.interfaceVar

    currentDensity = metalVar / metalConcentration * exchangeCurrentDensity

    depositionRateVariable = currentDensity * molarVolume / charge / faradaysConstant

    extensionVelocityVariable = CellVariable(
        name = 'extension velocity',
        mesh = mesh,
        value = depositionRateVariable)

    catalystSurfactantEquation = AdsorbingSurfactantEquation(
        catalystVar,
        distanceVar = distanceVar,
        bulkVar = 0,
        rateConstant = 0,
        consumptionCoeff = consumptionRateConstant * extensionVelocityVariable)

    advectionEquation = TransientTerm() + FirstOrderAdvectionTerm(extensionVelocityVariable)

    metalEquation = buildMetalIonDiffusionEquation(
        ionVar = metalVar,
        distanceVar = distanceVar,
        depositionRate = depositionRateVariable,
        diffusionCoeff = metalDiffusion,
        metalIonMolarVolume = molarVolume)

    metalVar.constrain(metalConcentration, mesh.facesTop)

    if displayViewers:

        try:
            from mayaviSurfactantViewer import MayaviSurfactantViewer
            viewer = MayaviSurfactantViewer(distanceVar, catalystVar.interfaceVar, zoomFactor = 1e6, datamax=1.0, datamin=0.0, smooth = 1, title = 'catalyst coverage', animate=True)

        except:

            class PlotVariable(CellVariable):
                def __init__(self, var = None, name = ''):
                    CellVariable.__init__(self, mesh = mesh.fineMesh, name = name)
                    self.var = self._requires(var)

                def _calcValue(self):
                    return numerix.array(self.var(self.mesh.cellCenters))

            viewer = MultiViewer(viewers=(
                    Viewer(PlotVariable(var = distanceVar), datamax=1e-9, datamin=-1e-9),
                    Viewer(PlotVariable(var = catalystVar.interfaceVar))))
    else:
        viewer = None

    levelSetUpdateFrequency = int(0.7 * narrowBandWidth / cellSize / cflNumber / 2)
    step = 0

    while step < numberOfSteps:

        if step % 10 == 0 and viewer is not None:
            viewer.plot()

        if step % levelSetUpdateFrequency == 0:

            distanceVar.calcDistanceFunction()

        extensionVelocityVariable.setValue(numerix.array(depositionRateVariable))

        dt = cflNumber * cellSize / max(extensionVelocityVariable.globalValue)
        distanceVar.extendVariable(extensionVelocityVariable)

        advectionEquation.solve(distanceVar, dt = dt)

        catalystSurfactantEquation.solve(catalystVar, dt = dt)

        metalEquation.solve(metalVar, dt = dt)

        step += 1

    point = ((5e-09,), (1.15e-07,))
    value = 1.45346701e-09
    return abs(float(distanceVar(point, order=1)) - value) < cellSize / 10.0

__all__ = ["runGold"]

if __name__ == '__main__':
    runGold(numberOfSteps = 300, cellSize = 0.05e-7, displayViewers=False)
    raw_input("finished")
