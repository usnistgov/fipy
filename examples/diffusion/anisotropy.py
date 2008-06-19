#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "circle.py"
 #                                    created: 4/6/06 {11:26:11 AM}
 #                                last update: 5/1/08 {3:25:37 PM}
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://ctcms.nist.gov
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2006- 4- 6 JEG 1.0 original
 # ###################################################################
 ##

r"""

This example demonstrates how to solve diffusion with an anisotropic
coefficient.  We wish to solve the following problem.

.. raw:: latex

    $$ \frac{\partial \phi}{\partial t} = \partial_j \Gamma_{ij}
    \partial_i \phi $$ on a circular domain centred at $(0, 0)$. We
    can choose an anisotropy of 80\% such that $$\Gamma' =
    \begin{pmatrix} 0.2 & 0 \\ 0 & 1 \end{pmatrix}$$ A new matrix is
    formed by rotating $\Gamma'$ such that $$R = \begin{pmatrix}
    \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta
    \end{pmatrix}$$ and $$ \Gamma = R \Gamma'
    R^T $$

    In the case of a point source at $(0, 0)$ a reference
    solution is given by, $$ \phi \left( X, Y, t \right) = Q \frac{
    \exp \left( -\frac{1}{4 t} \left( \frac{ X^2 }{ \Gamma'_{00}} +
    \frac{ Y^2 }{ \Gamma'_{11}} \right) \right) }{ 4 \pi t
    \sqrt{\Gamma'_{00} \Gamma'_{11}} } $$ where $ \left(X, Y \right)^T
    = R \left(x, y \right)^T $ and $Q$ is the initial mass.

..

    >>> from fipy import *

Import a mesh previously created using Gmsh.

    >>> mesh = GmshImporter2D(os.path.splitext(__file__)[0] + '.msh')

Set the center most cell to have a value.

    >>> var = CellVariable(mesh=mesh, hasOld=1)
    >>> x, y = mesh.getCellCenters()
    >>> var[numerix.argmin(x**2 + y**2)] = 1.

Choose an orientation for the anisotropy.

    >>> theta = numerix.pi / 4.
    >>> rotationMatrix = numerix.array(((numerix.cos(theta), numerix.sin(theta)), \
    ...                                 (-numerix.sin(theta), numerix.cos(theta))))
    >>> gamma_prime = numerix.array(((0.2, 0.), (0., 1.)))
    >>> DOT = numerix.NUMERIX.dot
    >>> gamma = DOT(DOT(rotationMatrix, gamma_prime), numerix.transpose(rotationMatrix))

Make the equation, viewer and solve.

    >>> eqn = TransientTerm() == DiffusionTerm((gamma,))

    >>> if __name__ == '__main__':
    ...     viewer = make(var, limits={'datamin' : 0.0, 'datamax' : 0.001})

    >>> mass = float(numerix.sum(mesh.getCellVolumes() * var))
    >>> time = 0
    >>> dt=0.00025 

    >>> for i in range(40):
    ...     var.updateOld()
    ...     res = 1.
    ...     
    ...     while res > 1e-2:
    ...         res = eqn.sweep(var, dt=dt)
    ...  
    ...     if __name__ == '__main__':
    ...         viewer.plot()
    ...     time += dt

Compare with the analytical solution (within 5% accuracy).

    >>> X, Y = numerix.dot(mesh.getCellCenters(), CellVariable(mesh=mesh, rank=2, value=rotationMatrix))
    >>> solution = mass * numerix.exp(-(X**2 / gamma_prime[0][0] + Y**2 / gamma_prime[1][1]) / (4 * time)) / (4 * numerix.pi * time * numerix.sqrt(gamma_prime[0][0] * gamma_prime[1][1]))
    >>> print max(abs((var - solution) / max(solution))) < 0.05
    True

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input('finished')
