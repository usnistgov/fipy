## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "__init__.py"
 #                                     created: 11/1/06 {4:27:40 PM}
 #                                 last update: 11/9/06 {7:44:28 PM}
 # Author: Jonathan Guyer
 # E-mail: <guyer@nist.gov>
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
 # History
 # 
 # modified   by  rev reason
 # ---------- --- --- -----------
 # 2006-11-06 JEG 1.0 original
 # 
 # ########################################################################
 ##

__docformat__ = 'restructuredtext'

from iterator import Iterator
from pseudoRKQSIterator import PseudoRKQSIterator
from pidIterator import PIDIterator

def L1norm(var, matrix, RHSvector):
    r"""
    :Parameters:
      - `var`: The `CellVariable` in question.
      - `matrix`: *(ignored)*
      - `RHSvector`: *(ignored)*
      
    :Returns: 
      .. raw:: latex

         \[
         \frac{\|\mathtt{var} - \mathtt{var}^\text{old}\|_1}{\|\mathtt{var}^\text{old}\|_1}
         \]
         where $\|\vec{x}\|_1 = \sum_{r=1}^{n} |x_r|$ is the $L^1$-norm of
         $\vec{x}$.
    """
    from fipy.tools.numerix import add
    denom = add.reduce(abs(var.getOld()))
    return add.reduce(abs(var - var.getOld())) / (denom + (denom == 0))
    
def L2norm(var, matrix, RHSvector):
    r"""
    :Parameters:
      - `var`: The `CellVariable` in question
      - `matrix`: *(ignored)*
      - `RHSvector`: *(ignored)*
      
    :Returns: 
      .. raw:: latex
    
         \[
         \frac{\|\mathtt{var} - \mathtt{var}^\text{old}\|_2}{\|\mathtt{var}^\text{old}\|_2}
         \]
         where $\|\vec{x}\|_2 = \sqrt{\sum_{r=1}^{n} |x_r|^2}$ is the
         $L^2$-norm of $\vec{x}$.
    """
    from fipy.tools.numerix import add, sqrt
    denom = sqrt(add.reduce(var.getOld()**2))
    return sqrt(add.reduce((var - var.getOld())**2)) / (denom + (denom == 0))
    
def sweepMonotonic(fn, *args, **kwargs):
    """
    Repeatedly calls `fn(*args, **kwargs)` until the residual returned by
    `fn()` is no longer decreasing.
    
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
