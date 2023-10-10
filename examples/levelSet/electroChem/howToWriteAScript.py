r"""Tutorial for writing an electrochemical superfill script.

This input file demonstrates how to
create a new superfill script if the existing suite of scripts do
not meet the required needs. It provides the functionality of
:mod:`examples.levelSet.electroChem.simpleTrenchSystem`.

To run this example from the base FiPy directory type::

    $ python examples/levelSet/electroChem/howToWriteAScript.py --numberOfElements=10000 --numberOfSteps=800

at the command line. The results of the simulation will be displayed
and the word ``finished`` in the terminal at the end of the
simulation. To obtain this example in a plain script file in order to
edit and run type::

    $ python setup.py copy_script --From examples/levelSet/electroChem/howToWriteAScript.py --To myScript.py

in the base :term:`FiPy` directory. The file :file:`myScript.py` will contain the
script.

The following is an explicit explanation of the input commands
required to set up and run the problem. At the top of the file all the
parameter values are set. Their use will be explained during the
instantiation of various objects and are the same as those explained in
:mod:`examples.levelSet.electroChem.simpleTrenchSystem`.

The following parameters (all in S.I. units)  represent,

 * physical constants,

   >>> faradaysConstant = 9.6e4
   >>> gasConstant = 8.314
   >>> transferCoefficient = 0.5

 * properties associated with the catalyst species,

   >>> rateConstant0 = 1.76
   >>> rateConstant3 = -245e-6
   >>> catalystDiffusion = 1e-9
   >>> siteDensity = 9.8e-6

 * properties of the cupric ions,

   >>> molarVolume = 7.1e-6
   >>> charge = 2
   >>> metalDiffusionCoefficient = 5.6e-10

 * parameters dependent on experimental constraints,

   >>> temperature = 298.
   >>> overpotential = -0.3
   >>> bulkMetalConcentration = 250.
   >>> catalystConcentration = 5e-3
   >>> catalystCoverage = 0.

 * parameters obtained from experiments on flat copper electrodes,

   >>> currentDensity0 = 0.26
   >>> currentDensity1 = 45.

 * general simulation control parameters,

   >>> cflNumber = 0.2
   >>> numberOfCellsInNarrowBand = 10
   >>> cellsBelowTrench = 10
   >>> cellSize = 0.1e-7

 * parameters required for a trench geometry,

   >>> trenchDepth = 0.5e-6
   >>> aspectRatio = 2.
   >>> trenchSpacing = 0.6e-6
   >>> boundaryLayerDepth = 0.3e-6

The hydrodynamic boundary layer depth (``boundaryLayerDepth``) is
intentionally small in this example to keep the mesh at a reasonable
size.

Build the mesh:

.. index::
   pair: module; fipy.tools.parser

>>> from fipy.tools.parser import parse
>>> numberOfElements = parse('--numberOfElements', action='store',
...     type='int', default=-1)
>>> numberOfSteps = parse('--numberOfSteps', action='store',
...     type='int', default=2)

.. index::
   single: sqrt
   single: exp

>>> from fipy import *

>>> if numberOfElements != -1:
...     pos = trenchSpacing * cellsBelowTrench / 4 / numberOfElements
...     sqr = trenchSpacing * (trenchDepth + boundaryLayerDepth) \
...           / (2 * numberOfElements)
...     cellSize = pos + numerix.sqrt(pos**2 + sqr)
... else:
...     cellSize = 0.1e-7

>>> yCells = cellsBelowTrench \
...          + int((trenchDepth + boundaryLayerDepth) / cellSize)
>>> xCells = int(trenchSpacing / 2 / cellSize)

.. index::
   single: Grid2D

>>> from .metalIonDiffusionEquation import buildMetalIonDiffusionEquation
>>> from .adsorbingSurfactantEquation import AdsorbingSurfactantEquation

>>> from fipy import serialComm
>>> mesh = Grid2D(dx=cellSize,
...               dy=cellSize,
...               nx=xCells,
...               ny=yCells,
...               communicator=serialComm)

A ``distanceVariable`` object,
:math:`\phi`, is  required to store  the  position of the interface.

The ``distanceVariable`` calculates its value so that it is a distance
function
(*i.e.* holds the distance at any point in the mesh from the electrolyte/metal
interface at :math:`\phi = 0`) and :math:`|\nabla\phi| = 1`.

First, create the :math:`\phi` variable, which is initially set to -1 everywhere.
Create an initial variable,

.. index::
   single: DistanceVariable

>>> narrowBandWidth = numberOfCellsInNarrowBand * cellSize
>>> distanceVar = DistanceVariable(
...    name='distance variable',
...    mesh= mesh,
...    value=-1.,
...    hasOld=1)

The electrolyte region will be the positive region of the domain while the metal
region will be negative.

>>> bottomHeight = cellsBelowTrench * cellSize
>>> trenchHeight = bottomHeight + trenchDepth
>>> trenchWidth = trenchDepth / aspectRatio
>>> sideWidth = (trenchSpacing - trenchWidth) / 2

>>> x, y = mesh.cellCenters
>>> distanceVar.setValue(1., where=(y > trenchHeight)
...                                 | ((y > bottomHeight)
...                                    & (x < xCells * cellSize - sideWidth)))

>>> distanceVar.calcDistanceFunction(order=2) #doctest: +LSM

The ``distanceVariable`` has now been created to mark the interface. Some other
variables need to be created that govern the concentrations of various species.

Create the catalyst surfactant coverage, :math:`\theta`, variable.
This variable influences the deposition rate.

.. index::
   single: SurfactantVariable

>>> catalystVar = SurfactantVariable(
...     name="catalyst variable",
...     value=catalystCoverage,
...     distanceVar=distanceVar)

Create the bulk catalyst concentration, :math:`c_{\theta}`,
in the electrolyte,

.. index::
   single: CellVariable

>>> bulkCatalystVar = CellVariable(
...     name='bulk catalyst variable',
...     mesh=mesh,
...     value=catalystConcentration)

Create the bulk metal ion concentration,
:math:`c_m`, in the electrolyte.

>>> metalVar = CellVariable(
...     name='metal variable',
...     mesh=mesh,
...     value=bulkMetalConcentration)

The following commands build the ``depositionRateVariable``,
:math:`v`. The ``depositionRateVariable`` is given by the following equation.

.. math::

   v = \frac{i \Omega}{n F}

where :math:`\Omega` is the metal molar volume, :math:`n` is the metal ion
charge and :math:`F` is Faraday's constant. The current density is given
by

.. math::

   i = i_0 \frac{c_m^i}{c_m^{\infty}} \exp{ \left( \frac{- \alpha F}{R T} \eta \right) }

where :math:`c_m^i` is the metal ion concentration in the bulk at the
interface, :math:`c_m^{\infty}` is the far-field bulk concentration of
metal ions, :math:`\alpha` is the transfer coefficient, :math:`R` is the gas
constant, :math:`T` is the temperature and :math:`\eta` is the
overpotential. The exchange current density is an empirical
function of catalyst coverage,

.. math::

   i_0(\theta) = b_0 + b_1 \theta

The commands needed to build this equation are,

>>> expoConstant = -transferCoefficient * faradaysConstant \
...                / (gasConstant * temperature)
>>> tmp = currentDensity1 \
...       * catalystVar.interfaceVar
>>> exchangeCurrentDensity = currentDensity0 + tmp
>>> expo = numerix.exp(expoConstant * overpotential)
>>> currentDensity = expo * exchangeCurrentDensity * metalVar \
...                  / bulkMetalConcentration
>>> depositionRateVariable = currentDensity * molarVolume \
...                          / (charge * faradaysConstant)

Build the extension velocity variable :math:`v_{\text{ext}}`. The extension
velocity uses the
``extensionEquation`` to spread the velocity at the interface to the
rest of the domain.

>>> extensionVelocityVariable = CellVariable(
...     name='extension velocity',
...     mesh=mesh,
...     value=depositionRateVariable)

Using the variables created above the governing equations will be
built.  The governing equation for surfactant conservation is given
by,

.. math::

   \dot{\theta} = J v \theta + k c_{\theta}^i (1 - \theta)

where :math:`\theta` is the coverage of catalyst at the interface,
:math:`J` is the curvature of the interface, :math:`v` is the normal velocity
of the interface, :math:`c_{\theta}^i` is the concentration of
catalyst in the bulk at the interface. The value :math:`k` is given
by an empirical function of overpotential,

.. math::

   k = k_0 + k_3 \eta^3

The above equation is represented by the
:class:`~examples.levelSet.electroChem.adsorbingSurfactantEquation.AdsorbingSurfactantEquation`
in :term:`FiPy`:

>>> surfactantEquation = AdsorbingSurfactantEquation(
...     surfactantVar=catalystVar,
...     distanceVar=distanceVar,
...     bulkVar=bulkCatalystVar,
...     rateConstant=rateConstant0 \
...                    + rateConstant3 * overpotential**3)

The variable :math:`\phi` is advected by the
``advectionEquation`` given by,

.. math::

   \frac{\partial \phi}{\partial t} + v_{\text{ext}}|\nabla \phi| = 0

and is set up with the following commands:

.. index::
   single: AdvectionTerm

>>> advectionEquation = TransientTerm() + AdvectionTerm(extensionVelocityVariable)

The diffusion of metal ions from the far field to the interface is
governed by,

.. math::

   \frac{\partial c_m}{\partial t} = \nabla \cdot D \nabla  c_m

where,

.. math::

   D = \begin{cases}
   D_m & \text{when $\phi > 0$,} \\
   0   & \text{when $\phi \le 0$}
   \end{cases}

The following boundary condition applies at :math:`\phi = 0`,

.. math::

   D \hat{n} \cdot \nabla c = \frac{v}{\Omega}.

The metal ion diffusion equation is set up with the following commands.

>>> metalEquation = buildMetalIonDiffusionEquation(
...     ionVar=metalVar,
...     distanceVar=distanceVar,
...     depositionRate=depositionRateVariable,
...     diffusionCoeff=metalDiffusionCoefficient,
...     metalIonMolarVolume=molarVolume,
... )

>>> metalVar.constrain(bulkMetalConcentration, mesh.facesTop)

The surfactant bulk diffusion equation solves the bulk diffusion of a
species with a source term for the jump from the bulk to an interface.
The governing equation is given by,

.. math::

   \frac{\partial c}{\partial t} = \nabla \cdot D \nabla  c

where,

.. math::

   D = \begin{cases}
   D_{\theta} & \text{when $\phi > 0$} \\
   0          & \text{when $\phi \le 0$}
   \end{cases}

The jump condition at the interface is defined by Langmuir
adsorption. Langmuir adsorption essentially states that the ability
for a species to jump from an electrolyte to an interface is
proportional to the concentration in the electrolyte, available site
density and a jump coefficient. The boundary condition
at :math:`\phi = 0` is given by,

.. math::

   D \hat{n} \cdot \nabla c = -k c (1 - \theta).

The surfactant bulk diffusion equation is set up with the following commands.

>>> from .surfactantBulkDiffusionEquation import buildSurfactantBulkDiffusionEquation
>>> bulkCatalystEquation = buildSurfactantBulkDiffusionEquation(
...     bulkVar=bulkCatalystVar,
...     distanceVar=distanceVar,
...     surfactantVar=catalystVar,
...     diffusionCoeff=catalystDiffusion,
...     rateConstant=rateConstant0 * siteDensity
... )

>>> bulkCatalystVar.constrain(catalystConcentration, mesh.facesTop)

If running interactively, create viewers.

.. index::
   single: MayaviSurfactantViewer

>>> if __name__ == '__main__':
...     try:
...         from .mayaviSurfactantViewer import MayaviSurfactantViewer
...         viewer = MayaviSurfactantViewer(distanceVar,
...                                         catalystVar.interfaceVar,
...                                         zoomFactor=1e6,
...                                         datamax=1.0,
...                                         datamin=0.0,
...                                         smooth=1)
...     except:
...         viewer = MultiViewer(viewers=(
...             Viewer(distanceVar, datamin=-1e-9, datamax=1e-9),
...             Viewer(catalystVar.interfaceVar)))
...         from fipy.models.levelSet.surfactant.matplotlibSurfactantViewer import MatplotlibSurfactantViewer
...         viewer = MatplotlibSurfactantViewer(catalystVar.interfaceVar)
... else:
...     viewer = None

The ``levelSetUpdateFrequency`` defines how often to call the
``distanceEquation`` to reinitialize the ``distanceVariable`` to a
distance function.

>>> levelSetUpdateFrequency = int(0.8 * narrowBandWidth \
...                               / (cellSize * cflNumber * 2))

The following loop runs for ``numberOfSteps`` time steps. The time step
is calculated with the CFL number and the maximum extension velocity.
:math:`v` to
:math:`v_\text{ext}` throughout the whole domain using
:math:`\nabla\phi\cdot\nabla v_\text{ext} = 0`.

>>> from builtins import range
>>> for step in range(numberOfSteps):
...
...     if viewer is not None:
...         viewer.plot()
...
...     if step % levelSetUpdateFrequency == 0:
...         distanceVar.calcDistanceFunction(order=2)
...
...     extensionVelocityVariable.setValue(depositionRateVariable())
...
...     distanceVar.updateOld()
...     distanceVar.extendVariable(extensionVelocityVariable, order=2)
...     dt = cflNumber * cellSize / extensionVelocityVariable.max()
...     advectionEquation.solve(distanceVar, dt=dt)
...     surfactantEquation.solve(catalystVar, dt=dt)
...     metalEquation.solve(var=metalVar, dt=dt)
...     bulkCatalystEquation.solve(var=bulkCatalystVar, dt=dt, solver=GeneralSolver()) #doctest: +LSM

The following is a short test case. It uses saved data from a
simulation with 5 time steps. It is not a test for accuracy but a way
to tell if something has changed or been broken.

.. index::
   single: loadtxt

>>> import os

>>> filepath = os.path.join(os.path.split(__file__)[0],
...                         "simpleTrenchSystem.gz")
>>> ##numerix.savetxt(filepath, numerix.array(catalystVar))
>>> print(catalystVar.allclose(numerix.loadtxt(filepath), rtol=1e-4)) #doctest: +LSMLIB
1

>>> from fipy import input
>>> if __name__ == '__main__':
...     input('finished')
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

def _run():
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript(__name__))

if __name__ == '__main__':
    _run()
