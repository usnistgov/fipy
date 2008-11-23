#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "howToWriteAScript.py"
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

r"""

.. raw:: latex

    \label{howToWriteAScript} This input file demonstrates how to
    create a new superfill script if the existing suite of scripts do
    not meet the required needs. It provides the functionality of
    Example~\ref{simpleTrenchSytem}.

To run this example from the base fipy directory type::
    
    $ examples/levelSet/electroChem/howToWriteAScript.py --numberOfElements=10000 --numberOfSteps=800

at the command line. The results of the simulation will be displayed
and the word `finished` in the terminal at the end of the
simulation. To obtain this example in a plain script file in order to
edit and run type::

    $ python setup.py copy_script --From examples/levelSet/electroChem/howToWriteAScript.py --To myScript.py

in the base FiPy directory. The file `myScript.py` will contain the
script.

The following is an explicit explanation of the input commands
required to set up and run the problem. At the top of the file all the
parameter values are set. Their use will be explained during the
instantiation of various objects and are the same as those explained in

.. raw:: latex

    Example~\ref{simpleTrenchSystem}.

The following parameters (all in S.I. units)  represent,

physical constants,

   >>> faradaysConstant = 9.6e4
   >>> gasConstant = 8.314
   >>> transferCoefficient = 0.5

properties associated with the catalyst species,

   >>> rateConstant0 = 1.76
   >>> rateConstant3 = -245e-6
   >>> catalystDiffusion = 1e-9
   >>> siteDensity = 9.8e-6
   
properties of the cupric ions,

   >>> molarVolume = 7.1e-6
   >>> charge = 2
   >>> metalDiffusionCoefficient = 5.6e-10

parameters dependent on experimental constraints,

   >>> temperature = 298.
   >>> overpotential = -0.3
   >>> bulkMetalConcentration = 250.
   >>> catalystConcentration = 5e-3
   >>> catalystCoverage = 0.
      
parameters obtained from experiments on flat copper electrodes,

   >>> currentDensity0 = 0.26
   >>> currentDensity1 = 45.

general simulation control parameters,

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

.. raw:: latex

   \IndexModule{parser}

..

   >>> from fipy.tools.parser import parse
   >>> numberOfElements = parse('--numberOfElements', action='store',
   ...     type='int', default=-1)
   >>> numberOfSteps = parse('--numberOfSteps', action='store',
   ...     type='int', default=5)

.. raw:: latex

   \IndexFunction{sqrt}
   \IndexFunction{exp}

..

   >>> if numberOfElements != -1:
   ...     pos = trenchSpacing * cellsBelowTrench / 4 / numberOfElements
   ...     sqr = trenchSpacing * (trenchDepth + boundaryLayerDepth) \
   ...           / (2 * numberOfElements)
   ...     cellSize = pos + sqrt(pos**2 + sqr)
   ... else:
   ...     cellSize = 0.1e-7

   >>> yCells = cellsBelowTrench \
   ...          + int((trenchDepth + boundaryLayerDepth) / cellSize)
   >>> xCells = int(trenchSpacing / 2 / cellSize)

.. raw:: latex

   \IndexClass{Grid2D}

..

   >>> from fipy import *
   >>> mesh = Grid2D(dx=cellSize,
   ...               dy=cellSize,
   ...               nx=xCells,
   ...               ny=yCells)

A `distanceVariable` object,

.. raw :: latex

    $\phi$, is  required to store  the  position of the interface  .

The `distanceVariable` calculates its value so that it is a distance
function 

.. raw:: latex

   (\emph{i.e.} holds the distance at any point in the mesh from the electrolyte/metal
   interface at $\phi$ = 0) and $|\nabla\phi| = 1$.

   First, create the $\phi$ variable, which is initially set to -1 everywhere. 
   Create an initial variable,
   \IndexClass{DistanceVariable}

..

   >>> narrowBandWidth = numberOfCellsInNarrowBand * cellSize
   >>> distanceVar = DistanceVariable(
   ...    name='distance variable',
   ...    mesh= mesh,
   ...    value=-1.,
   ...    narrowBandWidth=narrowBandWidth,
   ...    hasOld=1)

The electrolyte region will be the positive region of the domain while the metal
region will be negative.

   >>> bottomHeight = cellsBelowTrench * cellSize
   >>> trenchHeight = bottomHeight + trenchDepth
   >>> trenchWidth = trenchDepth / aspectRatio
   >>> sideWidth = (trenchSpacing - trenchWidth) / 2
   
   >>> x, y = mesh.getCellCenters()
   >>> distanceVar.setValue(1., where=(y > trenchHeight) 
   ...                                 | ((y > bottomHeight) 
   ...                                    & (x < xCells * cellSize - sideWidth)))

   >>> distanceVar.calcDistanceFunction(narrowBandWidth=1e10)

The `distanceVariable` has now been created to mark the interface. Some other
variables need to be created that govern the concentrations of various species.

.. raw:: latex

    Create the catalyst surfactant coverage, $\theta$, variable.
    This variable influences the deposition rate.
    \IndexClass{SurfactantVariable}

 ..

   >>> catalystVar = SurfactantVariable(
   ...     name="catalyst variable",
   ...     value=catalystCoverage,
   ...     distanceVar=distanceVar)

.. raw:: latex

    Create the bulk catalyst concentration, $c_{\theta}$,
    in the electrolyte,
    \IndexClass{CellVariable}

..

   >>> bulkCatalystVar = CellVariable(
   ...     name='bulk catalyst variable',
   ...     mesh=mesh,
   ...     value=catalystConcentration)
   
Create the bulk metal ion concentration,

.. raw:: latex

    $c_m$,

in the electrolyte.
        
   >>> metalVar = CellVariable(
   ...     name='metal variable',
   ...     mesh=mesh,
   ...     value=bulkMetalConcentration)

The following commands build the `depositionRateVariable`,

.. raw:: latex

    $v$.

The `depositionRateVariable` is given by the following equation.

.. raw:: latex

    $$ v = \frac{i \Omega}{n F} $$

    where $\Omega$ is the metal molar volume, $n$ is the metal ion
    charge and $F$ is Faraday's constant. The current density is given
    by

    $$ i = i_0 \frac{c_m^i}{c_m^{\infty}} \exp{ \left( \frac{- \alpha F}{R T} \eta \right) } $$

    where $c_m^i$ is the metal ion concentration in the bulk at the
    interface, $c_m^{\infty}$ is the far-field bulk concentration of
    metal ions, $\alpha$ is the transfer coefficient, $R$ is the gas
    constant, $T$ is the temperature and $\eta$ is the
    overpotential. The exchange current density is an empirical
    function of catalyst coverage,

    $$ i_0(\theta) = b_0 + b_1 \theta $$

The commands needed to build this equation are,

   >>> expoConstant = -transferCoefficient * faradaysConstant \
   ...                / (gasConstant * temperature)
   >>> tmp = currentDensity1 \
   ...       * catalystVar.getInterfaceVar()
   >>> exchangeCurrentDensity = currentDensity0 + tmp
   >>> expo = exp(expoConstant * overpotential)
   >>> currentDensity = expo * exchangeCurrentDensity * metalVar \
   ...                  / bulkMetalConcentration
   >>> depositionRateVariable = currentDensity * molarVolume \
   ...                          / (charge * faradaysConstant)

.. raw:: latex

    Build the extension velocity variable $v_{\text{ext}}$. The extension
    velocity uses the

`extensionEquation` to spread the velocity at the interface to the
rest of the domain.

   >>> extensionVelocityVariable = CellVariable(
   ...     name='extension velocity',
   ...     mesh=mesh,
   ...     value=depositionRateVariable)   

Using the variables created above the governing equations will be
built.  The governing equation for surfactant conservation is given
by,

.. raw:: latex

    $$ \dot{\theta} = J v \theta + k c_{\theta}^i (1 - \theta) $$

    where $\theta$ is the coverage of catalyst at the interface,
    $J$ is the curvature of the interface, $v$ is the normal velocity
    of the interface, $c_{\theta}^i$ is the concentration of
    catalyst in the bulk at the interface. The value $k$ is given
    by an empirical function of overpotential,
    $$ k = k_0 + k_3 \eta^3 $$
    The above equation is represented by the \Class{AdsorbingSurfactantEquation}
    in \FiPy{}:

..

   >>> surfactantEquation = AdsorbingSurfactantEquation(
   ...     surfactantVar=catalystVar,
   ...     distanceVar=distanceVar,
   ...     bulkVar=bulkCatalystVar,
   ...     rateConstant=rateConstant0 \
   ...                    + rateConstant3 * overpotential**3)

.. raw:: latex

    The variable $\phi$ is advected by the

`advectionEquation` given by,

.. raw:: latex

    $$ \frac{\partial \phi}{\partial t} + v_{\text{ext}}|\nabla \phi| = 0 $$
    and is set up with the following commands:
    \IndexFunction{buildHigherOrderAdvectionEquation}

..

   >>> advectionEquation = buildHigherOrderAdvectionEquation(
   ...     advectionCoeff=extensionVelocityVariable)

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
    The \verb|MetalIonDiffusionEquation| is set up with the following commands.
    \IndexClass{FixedValue}
    \IndexFunction{buildMetalIonDiffusionEquation}

..

   >>> metalEquation = buildMetalIonDiffusionEquation(
   ...     ionVar=metalVar,
   ...     distanceVar=distanceVar,
   ...     depositionRate=depositionRateVariable,
   ...     diffusionCoeff=metalDiffusionCoefficient,
   ...     metalIonMolarVolume=molarVolume,
   ... )

   >>> metalEquationBCs = FixedValue(faces=mesh.getFacesTop(), value=bulkMetalConcentration)

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
    The \verb|SurfactantBulkDiffusionEquation| is set up with the following commands.
    \IndexFunction{buildSurfactantBulkDiffusionEquation}

..

   >>> bulkCatalystEquation = buildSurfactantBulkDiffusionEquation(
   ...     bulkVar=bulkCatalystVar,
   ...     distanceVar=distanceVar,
   ...     surfactantVar=catalystVar,
   ...     diffusionCoeff=catalystDiffusion,
   ...     rateConstant=rateConstant0 * siteDensity
   ... )

   >>> catalystBCs = FixedValue(faces=mesh.getFacesTop(), value=catalystConcentration)
   
If running interactively, create viewers.

.. raw:: latex

   \IndexClass{MayaviSurfactantViewer}
   
..

   >>> if __name__ == '__main__':
   ...     try:
   ...         viewer = MayaviSurfactantViewer(distanceVar,
   ...                                         catalystVar.getInterfaceVar(),
   ...                                         zoomFactor=1e6,
   ...                                         datamax=1.0, 
   ...                                         datamin=0.0,
   ...                                         smooth=1)
   ...     except:
   ...         viewer = MultiViewer(viewers=(
   ...             Viewer(distanceVar, datamin=-1e-9, datamax=1e-9),
   ...             Viewer(catalystVar.getInterfaceVar())))
   ... else:
   ...     viewer = None

The `levelSetUpdateFrequency` defines how often to call the
`distanceEquation` to reinitialize the `distanceVariable` to a
distance function.

   >>> levelSetUpdateFrequency = int(0.8 * narrowBandWidth \
   ...                               / (cellSize * cflNumber * 2))

The following loop runs for `numberOfSteps` time steps. The time step
is calculated with the CFL number and the maximum extension velocity.

.. raw:: latex

    $v$ to
    $v_\text{ext}$ throughout the whole domain using
    $\nabla\phi\cdot\nabla v_\text{ext} = 0$.

..
   
   >>> for step in range(numberOfSteps):
   ...
   ...     if viewer is not None:
   ...         viewer.plot()
   ...
   ...     if step % levelSetUpdateFrequency == 0:
   ...         distanceVar.calcDistanceFunction()
   ...
   ...     extensionVelocityVariable.setValue(depositionRateVariable())
   ...
   ...     distanceVar.updateOld()
   ...     catalystVar.updateOld()
   ...     metalVar.updateOld()
   ...     bulkCatalystVar.updateOld()
   ...     distanceVar.extendVariable(extensionVelocityVariable)
   ...     dt = cflNumber * cellSize / extensionVelocityVariable.max()
   ...     advectionEquation.solve(distanceVar, dt=dt, solver=LinearCGSSolver())
   ...     surfactantEquation.solve(catalystVar, dt=dt)
   ...     metalEquation.solve(var=metalVar, dt=dt,
   ...                         boundaryConditions=metalEquationBCs, solver=LinearCGSSolver())
   ...     bulkCatalystEquation.solve(var=bulkCatalystVar, dt=dt,
   ...                                   boundaryConditions=catalystBCs, solver=LinearCGSSolver())

The following is a short test case. It uses saved data from a
simulation with 5 time steps. It is not a test for accuracy but a way
to tell if something has changed or been broken.

.. raw:: latex

   \IndexFunction{loadtxt}
   
..
 
   >>> import os
   >>> filepath = os.path.join(os.path.split(__file__)[0], 
   ...                         "simpleTrenchSystem.gz")
   >>> print catalystVar.allclose(loadtxt(filepath), rtol=1e-4)
   1

   >>> if __name__ == '__main__':
   ...     raw_input('finished')

"""
__docformat__ = 'restructuredtext'

def _run():
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript(__name__))
    
if __name__ == '__main__':
    _run()


