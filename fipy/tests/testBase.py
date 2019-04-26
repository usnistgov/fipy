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
