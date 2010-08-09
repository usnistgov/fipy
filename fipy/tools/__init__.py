#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "__init__.py"
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

try:
    import scipy
except:
    pass

try:
    from PyTrilinos import Epetra
    from fipy.tools.commWrapper import CommWrapper

    parallel = CommWrapper(Epetra=Epetra)

    if parallel.Nproc > 1:

        try:
            from mpi4py import MPI
            from fipy.tools.mpi4pyCommWrapper import Mpi4pyCommWrapper
            parallel = Mpi4pyCommWrapper(Epetra=Epetra, MPI=MPI)
        except ImportError:
            raise Exception("Could not import mpi4py. The package mpi4py is a required package if you are using Trilinos in parallel. Try installing using 'easy_install mpi4py'.")

    from fipy.tools.serialCommWrapper import SerialCommWrapper
    serial = SerialCommWrapper(Epetra=Epetra)

except ImportError:
    from fipy.tools.dummyComm import DummyComm
    parallel = DummyComm()
    serial = DummyComm()

import dump
import numerix
import vector
from dimensions.physicalField import PhysicalField
from numerix import *
from vitals import Vitals

