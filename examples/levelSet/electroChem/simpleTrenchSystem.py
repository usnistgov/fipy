#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "inputSimpleTrenchSystem.py"
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

r"""Model electrochemical superfill using the CEAC mechanism.

This input file
is a demonstration of the use of :term:`FiPy`
for modeling electrodeposition using the CEAC mechanism. The
material properties and experimental parameters used are roughly
those that have been previously
published :cite:`NIST:damascene:2003]`.

To run this example from the base fipy directory type::

    $ python examples/levelSet/electroChem/simpleTrenchSystem.py

at the command line. The results of the simulation will be displayed and the word
`finished` in the terminal at the end of the simulation. To run with a different
number of time steps change the ``numberOfSteps`` argument as follows,

.. index:: runSimpleTrenchSystem

>>> runSimpleTrenchSystem(numberOfSteps=2, displayViewers=False) #doctest: +LSMLIB
1

Change the ``displayViewers`` argument to ``True`` if you wish to see the
results displayed on the
screen. Example :mod:`examples.levelSet.electroChem.simpleTrenchSystem` gives explanation for
writing new scripts or modifying existing scripts that are
encapsulated by functions.

Any argument parameter can be changed. For example if the initial
catalyst coverage is not 0, then it can be reset,

>>> runSimpleTrenchSystem(numberOfSteps=2, catalystCoverage=0.1, displayViewers=False)#doctest: +LSM
0

The following image shows a schematic of a trench geometry along with
the governing equations for modeling electrodeposition with the CEAC
mechanism. All of the given equations are implemented in the
:func:`examples.levelSet.electroChem.simpleTrenchSystem.runSimpleTrenchSystem`
function. As stated above, all the parameters
in the equations can be changed with function arguments.

.. image:: electroChem/schematicOfEquations.*
   :width: 90%
   :align: center
   :alt: schematic of superfill equations

The following table shows the symbols used in the governing equations
and their corresponding arguments to the
:func:`~examples.levelSet.electroChem.simpleTrenchSystem.runSimpleTrenchSystem`
function. The boundary layer depth is intentionally small in this
example in order not to complicate the mesh. Further examples will
simulate more realistic boundary layer depths but will also have more
complex meshes requiring the :command:`gmsh` software.

.. index:: gmsh

.. this is kind of nasty, but reST tables can't handle what we need, particularly decimal alignment

.. math::

   \mbox{
    \begin{tabular}{|rllr@{.}ll|}
    \hline
    Symbol                & Description                       & Keyword Argument                      & \multicolumn{2}{l}{Value} & Unit                               \\
    \hline
    \multicolumn{6}{|c|}{Deposition Rate Parameters}                                                                                                                   \\
    \hline
    $v$                   & deposition rate                   &                                       & \multicolumn{2}{l}{}      & m s$^{-1}$                         \\
    $i$                   & current density                   &                                       & \multicolumn{2}{l}{}      & A m$^{-2}$                         \\
    $\Omega$              & molar volume                      & \texttt{molarVolume}                  & 7&1$\times$10$^{-6}$      & m$^3$ mol$^{-1}$                   \\
    $n$                   & ion charge                        & \texttt{charge}                       & \multicolumn{2}{c}{2}     &                                    \\
    $F$                   & Faraday's constant                & \texttt{faradaysConstant}             & 9&6$\times$10$^{-4}$      & C mol$^{-1}$                       \\
    $i_0$                 & exchange current density          &                                       & \multicolumn{2}{l}{}      & A m$^{-2}$                         \\
    $\alpha$              & transfer coefficient              & \texttt{transferCoefficient}          & 0&5                       &                                    \\
    $\eta$                & overpotential                     & \texttt{overpotential}                & -0&3                      & V                                  \\
    $R$                   & gas constant                      & \texttt{gasConstant}                  & 8&314                     & J K$^{-1}$ mol$^{-1}$               \\
    $T$                   & temperature                       & \texttt{temperature}                  & 298&0                     & K                                  \\
    $b_0$                 & current density for $\theta^0$    & \texttt{currentDensity0}              & 0&26                      & A m$^{-2}$                         \\
    $b_1$                 & current density for $\theta$      & \texttt{currentDensity1}              & 45&0                      & A m$^{-2}$                         \\
    \hline
    \multicolumn{6}{|c|}{Metal Ion Parameters}                                                                                                                         \\
    \hline
    $c_m$                 & metal ion concentration           & \texttt{metalConcentration}           & 250&0                     & mol m$^{-3}$                       \\
    $c_m^{\infty}$        & far field metal ion concentration & \texttt{metalConcentration}           & 250&0                     & mol m$^{-3}$                       \\
    $D_m$                 & metal ion diffusion coefficient   & \texttt{metalDiffusion}               & 5&6$\times$10$^{-10}$     & m$^2$ s$^{-1}$                     \\
    \hline
    \multicolumn{6}{|c|}{Catalyst Parameters}                                                                                                                          \\
    \hline
    $\theta$              & catalyst surfactant concentration & \texttt{catalystCoverage}             & 0&0                       &                                    \\
    $c_{\theta}$          & bulk catalyst concentration       & \texttt{catalystConcentration}        & 5&0$\times$10$^{-3}$      & mol m$^{-3}$                       \\
    $c_{\theta}^{\infty}$ & far field catalyst concentration  & \texttt{catalystConcentration}        & 5&0$\times$10$^{-3}$      & mol m$^{-3}$                       \\
    $D_{\theta}$          & catalyst diffusion coefficient    & \texttt{catalystDiffusion}            & 1&0$\times$10$^{-9}$      & m$^2$ s$^{-1}$                     \\
    $\Gamma$              & catalyst site density             & \texttt{siteDensity}                  & 9&8$\times$10$^{-6}$      & mol m$^{-2}$                       \\
    $k$                   & rate constant                     &                                       & \multicolumn{2}{l}{}      & m$^3$ mol$^{-1}$ s$^{-1}$          \\
    $k_0$                 & rate constant for $\eta^0$        & \texttt{rateConstant0}                & 1&76                      & m$^3$ mol$^{-1}$ s$^{-1}$          \\
    $k_3$                 & rate constant for $\eta^3$        & \texttt{rateConstant3}                & -245&0$\times$10$^{-6}$   & m$^3$ mol$^{-1}$ s$^{-1}$ V$^{-3}$ \\
    \hline
    \multicolumn{6}{|c|}{Geometry Parameters}                                                                                                                          \\
    \hline
    $D$                   & trench depth                      & \texttt{trenchDepth}                  & 0&5$\times$10$^{-6}$      & m                                  \\
    $D / W$               & trench aspect ratio               & \texttt{aspectRatio}                  & 2&0                       &                                    \\
    $S$                   & trench spacing                    & \texttt{trenchSpacing}                & 0&6$\times$10$^{-6}$      & m                                  \\
    $\delta$              & boundary layer depth              & \texttt{boundaryLayerDepth}           & 0&3$\times$10$^{-6}$      & m                                  \\
    \hline
    \multicolumn{6}{|c|}{Simulation Control Parameters}                                                                                                                \\
    \hline
                          & computational cell size           & \texttt{cellSize}                     & 0&1$\times$10$^{-7}$      & m                                  \\
                          & number of time steps              & \texttt{numberOfSteps}                & \multicolumn{2}{c}{5}     &                                    \\
                          & whether to display the viewers    & \texttt{displayViewers}               & \multicolumn{2}{c}{\texttt{True}} &                           \\
    \hline
    \end{tabular}
   }

If the MayaVi plotting
software is
installed (see :ref:`INSTALLATION`) then a plot should
appear that is updated every 20 time steps and will eventually
resemble the image below.

.. image:: electroChem/inputSimpleTrenchSystem.*
   :width: 90%
   :align: center
   :alt: catalyst coverage as a function of time during superfill

.. .. bibmissing:: /documentation/refs.bib
    :sort:
"""
__docformat__ = 'restructuredtext'

from fipy import CellVariable, DistanceVariable, SurfactantVariable, Grid2D, TransientTerm, AdvectionTerm, GeneralSolver, Viewer, MultiViewer
from fipy.tools import numerix
from metalIonDiffusionEquation import buildMetalIonDiffusionEquation
from adsorbingSurfactantEquation import AdsorbingSurfactantEquation

def runSimpleTrenchSystem(faradaysConstant=9.6e4,
                          gasConstant=8.314,
                          transferCoefficient=0.5,
                          rateConstant0=1.76,
                          rateConstant3=-245e-6,
                          catalystDiffusion=1e-9,
                          siteDensity=9.8e-6,
                          molarVolume=7.1e-6,
                          charge=2,
                          metalDiffusion=5.6e-10,
                          temperature=298.,
                          overpotential=-0.3,
                          metalConcentration=250.,
                          catalystConcentration=5e-3,
                          catalystCoverage=0.,
                          currentDensity0=0.26,
                          currentDensity1=45.,
                          cellSize=0.1e-7,
                          trenchDepth=0.5e-6,
                          aspectRatio=2.,
                          trenchSpacing=0.6e-6,
                          boundaryLayerDepth=0.3e-6,
                          numberOfSteps=5,
                          displayViewers=True):

    cflNumber = 0.2
    numberOfCellsInNarrowBand = 10
    cellsBelowTrench = 10

    yCells = cellsBelowTrench \
             + int((trenchDepth + boundaryLayerDepth) / cellSize)

    xCells = int(trenchSpacing / 2 / cellSize)

    from fipy.tools import serialComm
    mesh = Grid2D(dx = cellSize,
                  dy = cellSize,
                  nx = xCells,
                  ny = yCells,
                  communicator=serialComm)

    narrowBandWidth = numberOfCellsInNarrowBand * cellSize

    distanceVar = DistanceVariable(
        name = 'distance variable',
        mesh = mesh,
        value = -1.,
        hasOld = 1)

    bottomHeight = cellsBelowTrench * cellSize
    trenchHeight = bottomHeight + trenchDepth
    trenchWidth = trenchDepth / aspectRatio
    sideWidth = (trenchSpacing - trenchWidth) / 2

    x, y = mesh.cellCenters
    distanceVar.setValue(1., where=(y > trenchHeight) | ((y > bottomHeight) & (x < xCells * cellSize - sideWidth)))

    distanceVar.calcDistanceFunction(order=2)

    catalystVar = SurfactantVariable(
        name = "catalyst variable",
        value = catalystCoverage,
        distanceVar = distanceVar)

    bulkCatalystVar = CellVariable(
        name = 'bulk catalyst variable',
        mesh = mesh,
        value = catalystConcentration)

    metalVar = CellVariable(
        name = 'metal variable',
        mesh = mesh,
        value = metalConcentration)

    expoConstant = -transferCoefficient * faradaysConstant \
                   / (gasConstant * temperature)

    tmp = currentDensity1 * catalystVar.interfaceVar

    exchangeCurrentDensity = currentDensity0 + tmp

    expo = numerix.exp(expoConstant * overpotential)
    currentDensity = expo * exchangeCurrentDensity * metalVar \
                     / metalConcentration

    depositionRateVariable = currentDensity * molarVolume \
                             / (charge * faradaysConstant)

    extensionVelocityVariable = CellVariable(
        name = 'extension velocity',
        mesh = mesh,
        value = depositionRateVariable)

    surfactantEquation = AdsorbingSurfactantEquation(
        surfactantVar = catalystVar,
        distanceVar = distanceVar,
        bulkVar = bulkCatalystVar,
        rateConstant = rateConstant0 + rateConstant3 * overpotential**3)

    advectionEquation = TransientTerm() + AdvectionTerm(extensionVelocityVariable)

    metalEquation = buildMetalIonDiffusionEquation(
        ionVar = metalVar,
        distanceVar = distanceVar,
        depositionRate = depositionRateVariable,
        diffusionCoeff = metalDiffusion,
        metalIonMolarVolume = molarVolume,
    )

    metalVar.constrain(metalConcentration, mesh.facesTop)

    from surfactantBulkDiffusionEquation import buildSurfactantBulkDiffusionEquation
    bulkCatalystEquation = buildSurfactantBulkDiffusionEquation(
        bulkVar = bulkCatalystVar,
        distanceVar = distanceVar,
        surfactantVar = catalystVar,
        diffusionCoeff = catalystDiffusion,
        rateConstant = rateConstant0 * siteDensity
    )

    bulkCatalystVar.constrain(catalystConcentration, mesh.facesTop)

    if displayViewers:
        try:
            from mayaviSurfactantViewer import MayaviSurfactantViewer
            viewer = MayaviSurfactantViewer(distanceVar, catalystVar.interfaceVar, zoomFactor = 1e6, datamax=0.5, datamin=0.0, smooth = 1, title = 'catalyst coverage')
        except:
            viewer = MultiViewer(viewers=(
                Viewer(distanceVar, datamin=-1e-9, datamax=1e-9),
                Viewer(catalystVar.interfaceVar)))
    else:
        viewer = None

    levelSetUpdateFrequency = int(0.8 * narrowBandWidth \
                                  / (cellSize * cflNumber * 2))

    for step in range(numberOfSteps):

        if step>5 and step % 5 == 0 and viewer is not None:
            viewer.plot()

        if step % levelSetUpdateFrequency == 0:
            distanceVar.calcDistanceFunction(order=2)

        extensionVelocityVariable.setValue(depositionRateVariable())

        distanceVar.updateOld()

        distanceVar.extendVariable(extensionVelocityVariable, order=2)
        dt = cflNumber * cellSize / extensionVelocityVariable.max()

        advectionEquation.solve(distanceVar, dt = dt)
        surfactantEquation.solve(catalystVar, dt = dt)
        metalEquation.solve(metalVar, dt = dt)
        bulkCatalystEquation.solve(bulkCatalystVar, dt = dt, solver=GeneralSolver(tolerance=1e-15, iterations=2000))


    try:
        import os
        filepath = os.path.splitext(__file__)[0] + '.gz'
        print catalystVar.allclose(numerix.loadtxt(filepath), rtol = 1e-4)
    except:
        return 0

__all__ = ["runSimpleTrenchSystem"]

if __name__ == '__main__':
    runSimpleTrenchSystem(numberOfSteps = 800, cellSize = 0.05e-7)
    raw_input("finished")
