#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 8/26/04 {10:29:10 AM} 
 #                                last update: 3/7/05 {1:55:41 PM} { 1:23:41 PM}
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

r"""

.. raw:: latex

    This input file is a demonstration of the use of FiPy for modeling
    copper electroplating.  The material properties and experimental
    parameters used are roughly those that have been previously
    published~\cite{NIST:damascene:2003}.

To run this example from the base fipy directory type::
    
    $ examples/diffusion/steadyState/mesh1D/input.py

.. raw:: latex

    at the command line. The simulation took about 5 minutes on a
    computer with a 2GHz Athlon CPU. The results of the simulation

will be displayed and the word `finished` in the terminal at the
end of the simulation. The Gist or Gnuplot package is required to view the
results as the simulation is being executed (see the installation

.. raw:: latex

    guide in chapter~\ref{chap:Installation}).

The following is an explicit explanation of the input commands
required to set up and run the problem. At the top of the file all the
parameter values are set. Their use will be explained during the
instantiation of various objects.

The following parameters (all in S.I. units)  represent,

physical constants,

   >>> faradaysConstant = 9.6e4
   >>> gasConstant = 8.314
   >>> transferCoefficient = 0.5

properties associated with the accelerator species,

   >>> rateConstant = 1.76
   >>> overpotentialDependence = -245e-6
   >>> acceleratorDiffusionCoefficient = 1e-9
   >>> siteDensity = 9.8e-6
   
properties of the cupric ions,

   >>> atomicVolume = 7.1e-6,
   >>> charge = 2
   >>> metalDiffusionCoefficient = 5.6e-10

parameters dependent on experimental constraints,

   >>> temperature = 298.
   >>> overpotential = -0.3
   >>> bulkMetalConcentration = 250.
   >>> bulkAcceleratorConcentration = 5e-3
   >>> initialAcceleratorCoverage = 0.
      
parameters obtained from experiments on flat copper electrodes,

   >>> constantCurrentDensity = 0.26
   >>> acceleratorDependenceCurrentDensity = 45.

general simulation control parameters,

   >>> numberOfSteps = 300
   >>> cflNumber = 0.2
   >>> numberOfCellsInNarrowBand = 10
   >>> cellsBelowTrench = 10
   >>> cellSize = 0.1e-7
   
parameters required for a trench geometry,

   >>> trenchDepth = 0.5e-6
   >>> aspectRatio = 2.
   >>> trenchSpacing = 0.6e-6
   >>> boundaryLayerDepth = 0.3e-6
   
The hydrodynamic boundary layer depth (`boundaryLayerDepth`) is
intentionally small in this example to keep the mesh at a reasonable
size.

Build the mesh:

   >>> yCells = cellsBelowTrench + int((trenchDepth + boundaryLayerDepth) / cellSize)
   >>> xCells = int(trenchSpacing / 2 / cellSize)
   >>> from fipy.meshes.grid2D import Grid2D
   >>> mesh = Grid2D(dx = cellSize,
   ...               dy = cellSize,
   ...               nx = xCells,
   ...               ny = yCells)

A `distanceVariable` object,

.. raw :: latex

    $\phi$, is  required to store  the  position of the interface  (at
    $\phi$  =  0).

The `distanceVariable` calculates its value so that it is a a distance
function (i.e.  holds the distance at any point in the mesh from the
electrolyte/metal interface).

.. raw:: latex

    $|\nabla\phi| = 1$.

    Firstly, create the $\phi$

variable:

This is initially set to -1 everywhere. The electrolyte region will be
the positive region of the domain while the metal region will be
negative. Create a function for returning cells that lie in the electrolyte
region (positive region).

   >>> bottomHeight = cellsBelowTrench * cellSize
   >>> trenchHeight = bottomHeight + trenchDepth
   >>> trenchWidth = trenchDepth / aspectRatio
   >>> sideWidth = (trenchSpacing - trenchWidth) / 2
   >>> def electrolyteFunc(cell):
   ...     x,y = cell.getCenter()    
   ...     if y > trenchHeight:
   ...         return 1
   ...     elif y < bottomHeight:
   ...         return 0
   ...     elif x < sideWidth:
   ...         return 0
   ...     else:
   ...         return 1

Get the positive cells by passing the function,

   >>> electrolyteCells = mesh.getCells(electrolyteFunc)

Create an initial array,

   >>> import Numeric
   >>> values = -Numeric.ones(mesh.getNumberOfCells(), 'd')
   >>> for cell in electrolyteCells:
   ...     values[cell.getID()] = 1

   >>> narrowBandWidth = numberOfCellsInNarrowBand * cellSize
   >>> from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable        
   >>> distanceVar = DistanceVariable(
   ...    name = 'distance variable',
   ...    mesh = mesh,
   ...    value = values,
   ...    narrowBandWidth = narrowBandWidth)
   >>> distanceVar.calcDistanceFunction(narrowBandWidth = 1e10)

The `distanceVariable` has now been created to mark the interface. Some other
variables need to be created that govern the concentrations of various species.

.. raw:: latex

    Create the leveler and accelerator surfactant coverage, ($\theta_l
    and $\theta_a$), variables.

The accelerator variable influences the deposition rate while the
leveler variable removes accelerator from the interface.

   >>> from fipy.models.levelSet.surfactant.surfactantVariable import SurfactantVariable
   >>> levelerVar = SurfactantVariable(
   ...     name = "accelerator variable",
   ...     value = initialAcceleratorCoverage,
   ...     distanceVar = distanceVar)
   >>> acceleratorVar = SurfactantVariable(
   ...     name = "accelerator variable",
   ...     value = initialAcceleratorCoverage,
   ...     distanceVar = distanceVar)

.. raw:: latex

    Create the bulk accelerator concentration, $c_{\theta}$,

in the electrolyte,

   >>> from fipy.variables.cellVariable import CellVariable
   >>> bulkAcceleratorVar = CellVariable(
   ...     name = 'bulk accelerator variable',
   ...     mesh = mesh,
   ...     value = bulkAcceleratorConcentration)
   
Create the bulk metal ion concentration,

.. raw:: latex

    $c_m$,

in the electrolyte.
        
   >>> from fipy.variables.cellVariable import CellVariable
   >>> metalVar = CellVariable(
   ...     name = 'metal variable',
   ...     mesh = mesh,
   ...     value = bulkMetalConcentration)

The following commands build the `depositionRateVariable`,

.. raw:: latex

    $v$.

The `depositionRateVariable` is given by the following equation.

.. raw:: latex

    $$ v = \frac{i \Omega}{n F} $$

    where $\Omega$ is the metal atomic volume, $n$ is the metal ion
    charge and $F$ is Faraday's constant. The current density is given by

    $$ i = i_0 \frac{c_m^i}{c_m^{\infty}} \exp{ \left( \frac{- \alpha F}{R T} \eta \right) } $$

    where $c_m^i$ is the metal ion concentration in the bulk at the
    interface, $c_m^{\infty}$ is the far-field bulk concentration of
    metal ions, $\alpha$ is the transfer coefficient, $R$ is the gas
    constant, $T$ is the temperature and $\eta$ is the
    overpotential. The exchange current density is an empirical
    function of accelerator coverage,

    $$ i_0(\theta) = b_0 + b_1 \theta $$

The commands needed to build this equation are,

   >>> expoConstant = -transferCoefficient * faradaysConstant / gasConstant / temperature
   >>> tmp = acceleratorDependenceCurrentDensity * acceleratorVar.getInterfaceVar()
   >>> exchangeCurrentDensity = constantCurrentDensity + tmp
   >>> expo = Numeric.exp(expoConstant * overpotential)
   >>> currentDensity = exchangeCurrentDensity * metalVar / bulkMetalConcentration * expo
   >>> depositionRateVariable = currentDensity * atomicVolume / charge / faradaysConstant

.. raw:: latex

    Build the extension velocity variable $v_{\text{ext}}$. The extension
    velocity uses the

`extensionEquation` to spread the velocity at the interface to the
rest of the domain.

   >>> extensionVelocityVariable = CellVariable(
   ...     name = 'extension velocity',
   ...     mesh = mesh,
   ...     value = depositionRateVariable)   

Using the variables created above the governing equations will be
built.  The governing equation for surfactant conservation is given
by,

.. raw:: latex

    $$ \dot{\theta} = J v \theta + k c_{\theta}^i (1 - \theta) $$

    where $\theta$ is the coverage of accelerator at the interface,
    $J$ is the curvature of the interface, $v$ is the normal velocity
    of the interface, $c_{\theta}^i$ is the concentration of
    accelerator in the bulk at the interface. The value $k$ is given
    by an empirical function of overpotential,

    $$ k = k_0 + k_3 \eta^3 $$

The above equation is represented by the `AdsorbingSurfactantEquation`
in FiPy:

   >>> from fipy.models.levelSet.surfactant.adsorbingSurfactantEquation \
   ...             import AdsorbingSurfactantEquation
   
   >>> levelerSurfactantEquation = AdsorbingSurfactantEquation(
   ...     levelerVar,
   ...     distanceVar = distanceVar,
   ...     bulkVar = bulkAcceleratorConcentration,
   ...     rateConstant = (rateConstant + overpotentialDependence * overpotential**3))
   
   >>> acceleratorSurfactantEquation = AdsorbingSurfactantEquation(
   ...     acceleratorVar,
   ...     distanceVar = distanceVar,
   ...     bulkVar = bulkAcceleratorConcentration,
   ...     rateConstant = rateConstant + overpotentialDependence * overpotential**3,
   ...     otherVar = levelerVar,
   ...     otherBulkVar = bulkAcceleratorConcentration,
   ...     otherRateConstant = (rateConstant + overpotentialDependence * overpotential**3))
   

.. raw:: latex

    The variable $\phi$ is advected by the

`advectionEquation` given by,

.. raw:: latex

    $$ \frac{\partial \phi}{\partial t} + v_{\text{ext}}|\nabla \phi| = 0 $$

and is set up with the following commands:

   >>> from fipy.models.levelSet.advection.higherOrderAdvectionEquation \
   ...                import HigherOrderAdvectionEquation
   >>> advectionEquation = HigherOrderAdvectionEquation(
   ...     distanceVar,
   ...     advectionCoeff = extensionVelocityVariable)

The diffusion of metal ions from the far field to the interface is
governed by,

.. raw:: latex

    $$ \frac{\partial c_m}{\partial t} = \nabla \cdot D \nabla  c_m $$

    where,

    $$ D = \begin{cases}
    D_m & \text{when $\phi > 0$,} \\
    0   & \text{when $\phi \le 0$}
    \end{cases} $$

    The following boundary condition applies at $\phi = 0$,

    $$ D \hat{n} \cdot \nabla c = \frac{v}{\Omega}. $$

The `MetalIonDiffusionEquation` is set up with the following commands.

   >>> from fipy.boundaryConditions.fixedValue import FixedValue
   >>> from fipy.models.levelSet.electroChem.metalIonDiffusionEquation \
   ...                      import MetalIonDiffusionEquation
   >>> metalEquation = MetalIonDiffusionEquation(
   ...     metalVar,
   ...     distanceVar = distanceVar,
   ...     depositionRate = depositionRateVariable,
   ...     diffusionCoeff = metalDiffusionCoefficient,
   ...     metalIonAtomicVolume = atomicVolume,
   ...     boundaryConditions = (
   ...         FixedValue(
   ...             mesh.getFacesTop(),
   ...             bulkMetalConcentration
   ...         ),
   ...     )
   ... )

The `SurfactantBulkDiffusionEquation` solves the bulk diffusion of a
species with a source term for the jump from the bulk to an interface.
The governing equation is given by,

.. raw:: latex

    $$ \frac{\partial c}{\partial t} = \nabla \cdot D \nabla  c $$

where,

.. raw:: latex

    $$ D = \begin{cases}
    D_{\theta} & \text{when $\phi > 0$} \\
    0          & \text{when $\phi \le 0$}
    \end{cases} $$

The jump condition at the interface is defined by Langmuir
adsorption. Langmuir adsorption essentially states that the ability
for a species to jump from an electrolyte to an interface is
proportional to the concentration in the electrolyte, available site
density and a jump coefficient. The boundary condition

.. raw:: latex

    at $\phi = 0$ is given by,

    $$ D \hat{n} \cdot \nabla c = -k c (1 - \theta). $$

The `SurfactantBulkDiffusionEquation` is set up with the following commands.

   >>> from fipy.models.levelSet.surfactant.surfactantBulkDiffusionEquation \
   ...                 import SurfactantBulkDiffusionEquation
   >>> bulkAcceleratorEquation = SurfactantBulkDiffusionEquation(
   ...     bulkAcceleratorVar,
   ...     distanceVar = distanceVar,
   ...     surfactantVar = acceleratorVar,
   ...     diffusionCoeff = acceleratorDiffusionCoefficient,
   ...     rateConstant = rateConstant * siteDensity,
   ...         boundaryConditions = (
   ...         FixedValue(
   ...             mesh.getFacesTop(),
   ...             bulkAcceleratorConcentration
   ...         ),
   ...     ) 
   ... )

The equations are now given to an `Iterator` object in the order that
they will be solved.

   >>> from fipy.iterators.iterator import Iterator
   >>> iterator = Iterator((advectionEquation,
   ...                      levelerSurfactantEquation,
   ...                      acceleratorSurfactantEquation,
   ...                      metalEquation,
   ...                      bulkAcceleratorEquation))

The function below is constructed to encapsulate the creation of the
viewers.

   >>> def buildViewers():
   ...    resolution = 3
   ...    cells = yCells * 2**(resolution-1)
   ...    return (
   ...        fipy.viewers.make(
   ...            vars = distanceVar,
   ...            limits = {
   ...                'xmin': 0,
   ...                'xmax': cells,
   ...                'ymin': 0,
   ...                'xymax': cells,
   ...                'datamin': -1e-8, 
   ...                'datamax': 1e-8}),
   ...        fipy.viewers.make(
   ...            vars = acceleratorVar.getInterfaceVar(),
   ...            limits = {
   ...                'xmin': 0,
   ...                'xmax': cells,
   ...                'ymin': 0,
   ...                'xymax': cells}),
   ...        fipy.viewers.make(
   ...            vars = levelerVar.getInterfaceVar(),
   ...            limits = {
   ...                'xmin': 0,
   ...                'xmax': cells,
   ...                'ymin': 0,
   ...                'xymax': cells}))

The `levelSetUpdateFrequency` defines how often to call the
`distanceEquation` to reinitialize the `distanceVariable` to a
distance function.

   >>> levelSetUpdateFrequency = int(0.8 * narrowBandWidth / cellSize / cflNumber / 2)

The following loop runs for `numberOfSteps` time steps. The time step
is calculated with the CFL number and the maximum extension velocity.

.. raw:: latex

    $v$ to
    $v_\text{ext}$ throughout the whole domain using
    $\nabla\phi\cdot\nabla v_\text{ext} = 0$.

The following code steps through the solution.

   >>> if __name__ == '__main__':
   ...     viewers = buildViewers()
   ...     for step in range(numberOfSteps):
   ...         print 'step',step
   ...
   ...         if step % levelSetUpdateFrequency == 0:
   ...             distanceVar.calcDistanceFunction()
   ...
   ...         extensionVelocityVariable.setValue(Numeric.array(depositionRateVariable))
   ...
   ...         argmax = Numeric.argmax(extensionVelocityVariable)
   ...         dt = cflNumber * cellSize / extensionVelocityVariable[argmax]
   ...
   ...         distanceVar.extendVariable(extensionVelocityVariable)
   ...
   ...         iterator.timestep(dt = dt)
   ...         
   ...         for viewer in viewers:
   ...             viewer.plot()
   ...     raw_input('finished')

The following is a short test case. It uses saved data from a
simulation with 5 time steps. It is not a test for accuracy but a way
to tell if something has changed or been broken.

   >>> for i in range(5):
   ...     distanceVar.extendVariable(extensionVelocityVariable)
   ...     iterator.timestep(dt = 0.1)

   >>> import os
   >>> testFile = 'test.gz'
   >>> import examples.levelSet.electroChem
   >>> import gzip
   >>> filepath = os.path.join(examples.levelSet.electroChem.__path__[0], testFile)
   >>> filestream = gzip.open(filepath,'r')
   >>> import cPickle
   >>> testData = cPickle.load(filestream)
   >>> filestream.close()

   >>> Numeric.allclose(Numeric.array(acceleratorVar), testData)
   1
          
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    raw_input("finished")
