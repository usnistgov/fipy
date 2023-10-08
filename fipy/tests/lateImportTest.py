"""Classes to enable accumulating tests without importing them

Prevent failure to import one test from stopping execution of other tests.
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

import unittest

class _LateImportTestCase(unittest.TestCase):
    """:class:`~unittest.TestCase` that delays importing until test execution
    """

    def __str__(self):
        return "import %s" % self.moduleName

    def _post__init__(self, moduleName, suite):
        self.moduleName = moduleName
        self.suite = suite

    def setUp(self):
        """
        See documentation of `__import__` for why
        this ugly hack is necessary
        """
        module = __import__(self.moduleName)
        components = self.moduleName.split('.')
        for component in components[1:]:
            module = getattr(module, component)

        self.suite.addTest(self._getTestSuite(module = module))

    def _getTestSuite(self, module):
        return module._suite()

    def runTest(self):
        return

class _LateImportTestSuite(unittest.TestSuite):
    """:class:`~unittest.TestSuite` that delays importing until test execution
    """

    def __init__(self, testModuleNames = (), base = '__main__'):
        unittest.TestSuite.__init__(self)
        self._addTestModules(moduleNames = testModuleNames, base = base)

    def _addTestModules(self, moduleNames = (), base = '__main__'):
        for moduleName in moduleNames:
            self._addTestModule(moduleName = moduleName, base = base)

    def _addTestModule(self, moduleName, base = '__main__', testClass = _LateImportTestCase):
        if base == '__main__':
            base = []
        else:
            base = base.split('.')[:-1]

        test = testClass()
        test._post__init__(".".join(base + [moduleName]), self)
        self.addTest(test)
