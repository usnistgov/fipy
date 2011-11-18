## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "__init__.py"
 #
 # Author: Jonathan Guyer <guyer@nist.gov>
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
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
 # ########################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.steppers.stepper import Stepper
from fipy.steppers.pseudoRKQSStepper import PseudoRKQSStepper
from fipy.steppers.pidStepper import PIDStepper

__all__ = ["L1error", "L2error", "LINFerror", "sweepMonotonic"]

def residual(var, matrix, RHSvector):
    r"""
    Determines the residual for the current solution matrix and variable.
    
    :Parameters:
      - `var`: The `CellVariable` in question, *prior* to solution.
      - `matrix`: The coefficient matrix at this step/sweep
      - `RHSvector`: The 
      
    :Returns: 
      .. math::

         \|\mathsf{L}\vec{x} - \vec{b}\|_\infty
         
      where :math:`\|\vec{\xi}\|_\infty` is the :math:`L^\infty`-norm of :math:`\vec{\xi}`.
    """
    from fipy.tools.numerix import array, LINFnorm
    
    Lx = matrix * array(var)
    return LINFnorm(Lx - RHSvector)
    
def error(var, matrix, RHSvector, norm):
    r"""
    :Parameters:
      - `var`: The `CellVariable` in question.
      - `matrix`: *(ignored)*
      - `RHSvector`: *(ignored)*
      - `norm`: A function that will normalize its `array` argument and return 
        a single number
      
    :Returns: 
      .. math::

         \frac{\|\mathtt{var} - \mathtt{var}^\text{old}\|_?}
         {\|\mathtt{var}^\text{old}\|_?}
         
      where :math:`\|\vec{x}\|_?` is the normalization of :math:`\vec{x}` provided
      by :func:`~fipy.steppers.norm`.
    """
    from fipy.tools.numerix import L1norm
    denom = L1norm(var.old)
    return L1norm(var - var.old) / (denom + (denom == 0))

def L1error(var, matrix, RHSvector):
    r"""
    :Parameters:
      - `var`: The `CellVariable` in question.
      - `matrix`: *(ignored)*
      - `RHSvector`: *(ignored)*
      
    :Returns: 
      .. math::

         \frac{\|\mathtt{var} - \mathtt{var}^\text{old}\|_1}
         {\|\mathtt{var}^\text{old}\|_1}

      where :math:`\|\vec{x}\|_1` is the :math:`L^1`-norm of :math:`\vec{x}`.
    """
    from fipy.tools.numerix import L1norm
    return error(var, matrix, RHSvector, L1norm)
      
def L2error(var, matrix, RHSvector):
    r"""
    :Parameters:
      - `var`: The `CellVariable` in question.
      - `matrix`: *(ignored)*
      - `RHSvector`: *(ignored)*
      
    :Returns: 
      .. math::

         \frac{\|\mathtt{var} - \mathtt{var}^\text{old}\|_2}
         {\|\mathtt{var}^\text{old}\|_2}

      where :math:`\|\vec{x}\|_2` is the :math:`L^2`-norm of :math:`\vec{x}`.
    """
    from fipy.tools.numerix import L2norm
    return error(var, matrix, RHSvector, L2norm)
      
def LINFerror(var, matrix, RHSvector):
    r"""
    :Parameters:
      - `var`: The `CellVariable` in question.
      - `matrix`: *(ignored)*
      - `RHSvector`: *(ignored)*
      
    :Returns: 
      .. math::

         \frac{\|\mathtt{var} - \mathtt{var}^\text{old}\|_\infty}
         {\|\mathtt{var}^\text{old}\|_\infty}

      where :math:`\|\vec{x}\|_\infty` is the :math:`L^\infty`-norm of :math:`\vec{x}`.
    """
    from fipy.tools.numerix import LINFnorm
    return error(var, matrix, RHSvector, LINFnorm)
             
def sweepMonotonic(fn, *args, **kwargs):
    """
    Repeatedly calls :func:`fn(*args, **kwargs)` until the residual returned by
    :func:`fn` is no longer decreasing.
    
    :Parameters:
      - `fn`: The function to call
      - `args`: The unnamed function argument `list`
      - `kwargs`: The named function argument `dict`
      
    :Returns: the final residual
    """
    oldres = 1e20
    # res is 0 the first time
    fn(*args, **kwargs)
    res = fn(*args, **kwargs)
    while res < oldres:
        oldres = res
        res = fn(*args, **kwargs)
        
    return res
