#!/usr/bin/env python

## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "benchmarker.py"
 #
 # Author: Jonathan Guyer <guyer@nist.gov>
 # Author: Daniel Wheeler <daniel.wheeler@nist.gov>
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
 # ########################################################################
 ##

import time

from fipy.tools.parser import parse

class Benchmarker:
    def __init__(self):
        self.measureMemory = parse('--measureMemory', action = 'store_true', default=False)
        self.sampleTime = parse('--sampleTime', action = 'store', type = 'float', default = 1)
        
        self.keys = ['mesh', 'variables', 'terms', 'solver', 'BCs', 'solve']
        
        self.times = {}
        for key in self.keys:
            self.times[key] = 0

        if self.measureMemory:
            from fipy.tools.memoryLogger import MemoryLogger
            self.logger = MemoryLogger(sampleTime = self.sampleTime)
            self.memories = {}
            self.logger.start()
##             time.sleep(3)
            self.memories['baseline'] = self.logger.stop()
            for key in self.keys:
                self.memories[key] = self.memories['baseline']

        self.t0 = time.clock()
        
        
    def start(self):
        if self.measureMemory:
            self.logger.start()
            
        self.t1 = time.clock()

    def stop(self, name):
        self.times[name] = time.clock() - self.t1
        if self.measureMemory:
            self.memories[name] = self.logger.stop()
            
    def report(self, numberOfElements=1, steps=1):
        self.times['total'] = time.clock() - self.t0

        output = []
        maxMemory = -1
        for key in self.keys:
            if self.measureMemory:
                memory = self.memories[key] - self.memories['baseline']
                output += [str(memory)]
                maxMemory = max(maxMemory, memory)
            else:
                output += [str(self.times[key])]
                
        if self.measureMemory:
            output += [str(maxMemory), str(float(maxMemory) / numberOfElements)]
        else:
            output += [str(self.times['total']), str(self.times['solve'] / steps / numberOfElements)]

        return "\t".join(output)