r"""Solve the diffusion equation with an anisotropic diffusion coefficient.

We wish to solve the problem

.. math::

   \frac{\partial \phi}{\partial t} = \partial_j \Gamma_{ij} \partial_i \phi

on a circular domain centered at :math:`(0, 0)`. We can choose an anisotropy ratio of 5
such that

.. math::

   \Gamma' = \begin{pmatrix}
       0.2 & 0 \\
         0 & 1
   \end{pmatrix}

A new matrix is formed by rotating :math:`\Gamma'` such that

.. math::

   R = \begin{pmatrix}
        \cos\theta & \sin\theta \\
       -\sin\theta & \cos\theta
   \end{pmatrix}

and

.. math::

   \Gamma = R \Gamma' R^T

In the case of a point source at :math:`(0, 0)` a reference
solution is given by,

.. math::

   \phi \left( X, Y, t \right) = Q \frac{
    \exp \left( -\frac{1}{4 t} \left( \frac{ X^2 }{ \Gamma'_{00}} +
    \frac{ Y^2 }{ \Gamma'_{11}} \right) \right) }{ 4 \pi t
    \sqrt{\Gamma'_{00} \Gamma'_{11}} }

where :math:`\left(X, Y \right)^T = R \left(x, y \right)^T` and :math:`Q` is the initial
mass.

>>> from fipy import CellVariable, Gmsh2D, Viewer, TransientTerm, DiffusionTermCorrection
>>> from fipy.tools import serialComm, numerix

Import a mesh previously created using :term:`Gmsh`.

>>> import os
>>> mesh = Gmsh2D(os.path.splitext(__file__)[0] + '.msh', communicator=serialComm) # doctest: +GMSH

Set the center-most cell to have a value.

>>> var = CellVariable(mesh=mesh, hasOld=1) # doctest: +GMSH
>>> x, y = mesh.cellCenters # doctest: +GMSH
>>> var[numerix.argmin(x**2 + y**2)] = 1. # doctest: +GMSH

Choose an orientation for the anisotropy.

>>> theta = numerix.pi / 4.
>>> rotationMatrix = numerix.array(((numerix.cos(theta), numerix.sin(theta)), \
...                                 (-numerix.sin(theta), numerix.cos(theta))))
>>> gamma_prime = numerix.array(((0.2, 0.), (0., 1.)))
>>> DOT = numerix.NUMERIX.dot
>>> gamma = DOT(DOT(rotationMatrix, gamma_prime), numerix.transpose(rotationMatrix))

Make the equation, viewer and solve.

>>> eqn = TransientTerm() == DiffusionTermCorrection((gamma,))

>>> if __name__ == '__main__':
...     viewer = Viewer(var, datamin=0.0, datamax=0.001)

>>> mass = float(var.cellVolumeAverage * numerix.sum(mesh.cellVolumes)) # doctest: +GMSH
>>> time = 0
>>> dt=0.00025

>>> from builtins import range
>>> for i in range(20):
...     var.updateOld() # doctest: +GMSH
...     res = 1.
... 
...     while res > 1e-2:
...         res = eqn.sweep(var, dt=dt) # doctest: +GMSH
... 
...     if __name__ == '__main__':
...         viewer.plot()
...     time += dt

Compare with the analytical solution (within 5% accuracy).

>>> X, Y = numerix.dot(mesh.cellCenters, CellVariable(mesh=mesh, rank=2, value=rotationMatrix)) # doctest: +GMSH
>>> solution = mass * numerix.exp(-(X**2 / gamma_prime[0][0] + Y**2 / gamma_prime[1][1]) / (4 * time)) / (4 * numerix.pi * time * numerix.sqrt(gamma_prime[0][0] * gamma_prime[1][1])) # doctest: +GMSH
>>> print(max(abs((var - solution) / max(solution))) < 0.08) # doctest: +GMSH
True

"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')
