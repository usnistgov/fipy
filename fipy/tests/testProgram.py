from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

import unittest

class TestProgram(unittest.TestProgram):
    """A command-line program that runs a set of tests

    This is primarily for making test modules conveniently executable.
    """
    def parseArgs(self, argv):
        import getopt
##      inline = 0
##        numMesh = 0
        try:
            options, args = getopt.getopt(argv[1:], 'hHvq',
                                          ['help', 'verbose', 'quiet', 'inline', 'Trilinos', 'Pysparse', 'pysparse', 'trilinos', 'no-pysparse', 'scipy', 'pyamg', 'skfmm', 'lsmlib'])
            for opt, value in options:
                if opt in ('-h', '-H', '--help'):
                    self.usageExit()
                if opt in ('-q', '--quiet'):
                    self.verbosity = 0
                if opt in ('-v', '--verbose'):
                    self.verbosity = 2
##              if opt in ('--inline',):
##                  inline = 1
##                if opt in ('--numMesh',):
##                    numMesh = 1
            if len(args) == 0 and self.defaultTest is None:
                self.test = self.testLoader.loadTestsFromModule(self.module)
                return
            if len(args) > 0:
                self.testNames = args
            else:
                self.testNames = (self.defaultTest,)
            self.createTests()
##            print argv
##            raw_input()
##          if inline:
##              argv[1:] = ['--inline']
##            if numMesh:
##                argv[1:] = ['--numMesh']
        except getopt.error as msg:
            self.usageExit(msg)

    def runTests(self):
        from fipy.tools import numerix
        printoptions = numerix.get_printoptions()
        if "legacy" in printoptions:
            numerix.set_printoptions(legacy="1.13")

        super(TestProgram, self).runTests()

main = TestProgram
