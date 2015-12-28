#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "efficiency_test.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Andrew Acquaviva <andrewa@nist.gov>
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

from distutils.core import Command
import os
import sys

__all__ = ["Efficiency_test"]

class Efficiency_test(Command):
    description = "run FiPy efficiency tests"

    user_options = [ ('minimumelements=', None, 'minimum number of elements'),
                     ('factor=', None, 'factor by which the number of elements is increased'),
                     ('inline', None, 'turn on inlining for the efficiency tests'),
                     ('cache', None, 'turn on variable caching'),
                     ('maximumelements=', None, 'maximum number of elements'),
                     ('sampleTime=', None, 'sampling interval for memory high-water'),
                     ('path=', None, 'directory to place output results in'),
                     ('uploadToCodespeed', None, 'flag to upload data to Codespeed'),
                     ('otherExample=', None, 'designate examples other than the default ones to benchmark'),
                     ('newElements=', None, 'alter the number of elements. NB: will not work on all examples'),
                     ('revisionNumber=', None, 'generate history of example(s) from given revision number')]

    def initialize_options(self):
        self.factor = 10
        self.inline = 0
        self.cache = 0
        self.maximumelements = 10000
        self.minimumelements = 100
        self.sampleTime = 1
        self.path = None
        self.otherExample = None
        self.newElements = None
        self.cases = ['examples/cahnHilliard/mesh2D.py', 'examples/phase/anisotropy.py', \
                          'examples/reactiveWetting/liquidVapor2D.py']
        self.uploadToCodespeed = False
        self.revisionNumber = None

    def finalize_options(self):
        self.factor = int(self.factor)
        self.maximumelements = int(self.maximumelements)
        self.minimumelements = int(self.minimumelements)
        self.sampleTime = float(self.sampleTime)
        if self.otherExample is not None:
            self.cases = [self.otherExample]

    def run(self):
        import time
        from fipy.tools.performance import efficiencyTestGenerator

        if self.revisionNumber is not None:
            self.revisionNumber = int(self.revisionNumber)
            from fipy.tools import efficiencyTestHistory
            import sys
            efficiencyTestHistory.run(self.revisionNumber)
            sys.exit()
        newCases = efficiencyTestGenerator.run(self.cases,self.newElements)

        for case in newCases:
            print "case: %s" % case

            if self.path is None:
                testPath = os.path.split(case)[0]
            else:
                testPath = self.path

            if not os.access(testPath, os.F_OK):
                os.makedirs(testPath)

            testPath = os.path.join(testPath, '%s.dat' % os.path.split(case)[1])

            if not os.path.isfile(testPath):
                f = open(testPath, 'w')

                f.write("\t".join(["--inline".center(10), "--cache".center(10), "Date".center(25),\
                        "Elements".center(10), "total runtime(s)".rjust(15)]))
                f.write("\n")
                f.flush()
            else:
                f = open(testPath, 'a')

            numberOfElements = self.minimumelements

            print "Running example:"
            cmd = ["python", "-W ignore", case]

            output = "\t".join([str(self.inline).center(10), str(self.cache).center(10),\
                                    (time.ctime()).center(25), str(numberOfElements).center(10)])
            print 'cmd',cmd
            w, r = os.popen2(cmd)

            outputlist= r.read().split()
            print outputlist
            init_time = outputlist[outputlist.index('Initialization-time:')+1]
            frst_timestp = outputlist[outputlist.index('First-timestep:')+1]
            avg_timestp = outputlist[outputlist.index('Average-timestep:')+1]
            runtime = outputlist[outputlist.index('Runtime:')+1]
            print "runtime: ", runtime
            output += '\t' + ''.join(runtime).strip()
            r.close()
            w.close()

            if numberOfElements == self.maximumelements:
                f.write(output + '\n' + "-"*100 + '\n')
                f.flush()
            else:
                f.write(output + '\n')
                f.flush()

            if self.uploadToCodespeed:
                import urllib, urllib2
                import time
                import pysvn
                from datetime import datetime

                CODESPEED_URL = 'http://build.cmi.kent.edu/codespeed/'

                revnum = pysvn.Client().info('.')['revision'].number
                revdate  = pysvn.Client().info('.')['commit_time']
                print 'revdate: ', datetime.fromtimestamp(revdate)

                def add(data):
                    params = urllib.urlencode(data)
                    response = "None"
                    print "Saving result for executable %s, revision %s, benchmark %s" % (
                        data['executable'], data['commitid'], data['benchmark'])
                    try:
                        f = urllib2.urlopen(CODESPEED_URL + 'result/add/', params)
                    except urllib2.HTTPError, e:
                        print str(e)
                        print e.read()
                        return
                    response = f.read()
                    f.close()
                    print "Server (%s) response: %s\n" % (CODESPEED_URL, response)

                benchmarks = ['Initialization', 'First timestep',\
                                  'Average of remaining timesteps', 'Total Runtime']
                results = [init_time, frst_timestp, avg_timestp, runtime]
                for i in range(len(benchmarks)):
                    data = {
                        'commitid': revnum,
                        'branch': 'default', #Always use default for trunk/master/tip
                        'project': 'FiPy',
                        'revision_date': datetime.fromtimestamp(revdate),
                        'executable': case,
                        'benchmark': benchmarks[i],
                        'environment': 'Sandbox A201590', ##"Debian A203166",
                        'result_value': results[i],
                        'result_date': datetime.fromtimestamp(revdate)
                        }
                    print datetime.fromtimestamp(revdate)
                    add(data)
            numberOfElements *= self.factor
            f.close()
            os.remove(case)
            os.remove("%s.dat" % case)
            print "Deleted temporary file ", case
