from distutils.core import Command
import os
import string
import sys
import glob
import ez_setup
from setuptools.command.test import test as _test

class Efficiency_test(Command):
    description = "run FiPy efficiency tests"

    user_options = [ ('minimumelements=', None, 'minimum number of elements'),
                     ('factor=', None, 'factor by which the number of elements is increased'),
                     ('inline', None, 'turn on inlining for the efficiency tests'),
                     ('cache', None, 'turn on variable caching'),
                     ('maximumelements=', None, 'maximum number of elements'),
                     ('sampleTime=', None, 'sampling interval for memory high-water'),
                     ('path=', None, 'directory to place output results in'),
                     ('uploadToCodespeed', None, 'flag to upload data to Codespeed')]
#                     ('example=', None, 'prompts the example to be benchmarked')]
    
    def initialize_options(self):
        w,r = os.popen2('python fipy/tools/generator.py')
        self.factor = 10
        self.inline = 0
        self.cache = 0
        self.maximumelements = 10000
        self.minimumelements = 100
        self.sampleTime = 1
        self.path = None
##        self.cases = ['examples/cahnHilliard/mesh2D.py']
        self.cases = ['./mesh1D.py']
        print self.cases
        self.uploadToCodespeed = False
#        self.example = "examples/cahnHilliard/mesh2D.py"
    
    def finalize_options(self):
        self.factor = int(self.factor)
        self.maximumelements = int(self.maximumelements)
        self.minimumelements = int(self.minimumelements)
        self.sampleTime = float(self.sampleTime)

    def run(self):
        import time
        import os

        for case in self.cases:
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

            while numberOfElements <= self.maximumelements:
                print "\tnumberOfElements: %i" % numberOfElements
                
##                cmd = ["python", "-W ignore", case, '--numberOfElements=%i' % numberOfElements, '--no-display nodisp']
                cmd = ["python", "-W ignore", case]

                output = "\t".join([str(self.inline).center(10), str(self.cache).center(10),\
                                   (time.ctime()).center(25), str(numberOfElements).center(10)])
                w, r = os.popen2(cmd)
                
##                timeCmd = cmd + ['--measureTime runtime']
##                w, r = os.popen2(' '.join(timeCmd))
##                print "' '.join(timeCmd): ", ' '.join(timeCmd)
##                raw_input()
                outputlist= r.read().split()
                runtime = outputlist[outputlist.index('runtime:')+1]
                print "runtime: ", runtime
                output += '\t' + ''.join(runtime).strip()
                r.close()
                w.close()

#                memCmd = cmd + ['--measureMemory', '--sampleTime=%f' % self.sampleTime]

#                w, r = os.popen4(' '.join(memCmd))
#                output += '\t' + ''.join(r.readlines()).strip()
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
      
                    CODESPEED_URL = "http://localhost:8000/"
                    revnum = pysvn.Client().info('.')['revision'].number
                    revdate  = pysvn.Client().info('.')['commit_time']

                    def add(data):
                        params = urllib.urlencode(data)
                        response = "None"
                        print "Executable %s, revision %s, benchmark %s" % (data['executable'],\
                                                            data['commitid'], data['benchmark']) 
                        g = urllib2.urlopen(CODESPEED_URL + 'result/add/', params)  
                        response = g.read()
                        g.close()
                        print "Server (%s) response: %s\n" % (CODESPEED_URL, response)
 
                    data = {
                        'commitid':revnum,
                        'branch': 'efficiency_test',#Always use default for trunk/master/tip
                        'project': 'FiPy',
                        'revision_date': datetime.fromtimestamp(revdate),
                        'executable': case,
                        'benchmark': 'float',
                        'environment': "FiPy",
                        'result_value': runtime,
                        'result_date': datetime.today(),
                        }                                     
                    add(data)   
                numberOfElements *= self.factor
                f.close
            
