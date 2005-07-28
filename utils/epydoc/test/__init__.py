# Natural Language Toolkit: Unit Tests
#
# Copyright (C) 2001 University of Pennsylvania
# Author: Edward Loper <edloper@gradient.cis.upenn.edu>
# URL: <http://nltk.sf.net>
# For license information, see LICENSE.TXT
#
# $Id$

"""
Unit tests for the NLTK modules.
"""

import unittest, sys, os, os.path

def testsuite():
    """
    Return a PyUnit testsuite for all tests.
    """
    test_modules = []
    
    import epydoc.test
    path = os.path.split(epydoc.test.__file__)[0]
    for f in os.listdir(path):
        if f[-3:] != '.py': continue
        if f[:9] == '__init__.': continue
        name = f[:-3]
        try:
            print 'Importing %s' % name
            exec('import %s as module' % name)
            test_modules.append(module)
        except KeyboardInterrupt: raise 
        except Exception, e: print 'Error importing %r: %s' % (name, e)
        except SystemExit, e: print 'Error importing %r: %s' % (name, e)
        except: print 'Error importing %r'

    return unittest.TestSuite([m.testsuite() for m in test_modules]) 

def test():
    """
    Run unit tests for the NLP toolkit; print results to stdout/stderr
    """

def usage(name):
    print """
    Usage: %s
    """ % name
    return 0

if __name__ == '__main__':
    if len(sys.argv) != 1: sys.exit(usage(sys.argv[0]))
    runner = unittest.TextTestRunner()
    runner.run(testsuite())

