#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
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
 # ###################################################################
 ##

r"""
How to update scripts from version 2.0 to 3.0.

:term:`FiPy` 3.0 introduces several syntax changes from :term:`FiPy` 2.0. We appreciate that
this is very inconvenient for our users, but we hope you'll agree that the new
syntax is easier to read and easier to use. We assure you that this is not
something we do casually; it has been over two and a half years since our last
incompatible change (when :term:`FiPy` 2.0 superceded :term:`FiPy` 1.0).

All examples included with version 3.0 have been updated to use the new syntax,
but any scripts you have written for :term:`FiPy` 2.0 will need to be updated. A
complete listing of the changes needed to take the :term:`FiPy` examples scripts from
version 2.0 to version 3.0 can be found at
    
    http://www.matforge.org/fipy/wiki/upgrade2_0examplesTo3_0
    
but we summarize the necessary changes here. If these tips are not sufficient to
make your scripts compatible with :term:`FiPy` 3.0, please don't hesitate to ask for
help on the `mailing list`_.


The following items **must** be changed in your scripts

 * We have reconsidered the change in FiPy 2.0 that included all of the 
   functions of the :mod:`~fipy.tools.numerix` module in the :mod:`fipy` namespace. 
   You now must be more explicit when referring to any of these functions:
       
   >>> from fipy import *
   >>> y = numerix.exp(x)

   >>> from fipy.tools.numerix import exp
   >>> y = exp(x)
       
   We generally use the first, but you may see us import specific functions if
   we feel it improves readability. You should feel free to use whichever form
   you find most comfortable.
   
   .. note:: the old behavior can be obtained, at least for now, but setting the
             :envvar:`FIPY_INCLUDE_NUMERIX_ALL` environment variable.
   
.. _mailing list:         http://www.ctcms.nist.gov/fipy/mail.html
"""
__docformat__ = 'restructuredtext'

