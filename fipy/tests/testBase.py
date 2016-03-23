#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "testBase.py"
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

__docformat__ = 'restructuredtext'

__all__ = []

import unittest

from fipy.tools import numerix

class _TestBase(unittest.TestCase):
    def assertWithinTolerance(self, first, second, tol = 1e-10, msg=None):
	"""Fail if the two objects are unequal by more than tol.
	"""
	if abs(first - second) > tol:
	    raise self.failureException(msg or '%s !~ %s' % (first, second))

    def assertArrayWithinTolerance(self, first, second, atol = 1e-10, rtol = 1e-10, msg=None):
	"""Fail if the two objects are unequal by more than tol.
	"""
	if not numerix.allclose(first, second, rtol = rtol, atol = atol):
	    raise self.failureException(msg or '\n%s\nis not\n%s' % (first, second))

    def testResult(self):
        pass
