## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "stepper.py"
 #                                     created: 10/31/06 {9:50:24 AM}
 #                                 last update: 11/10/06 {4:21:51 PM}
 # Author: Jonathan Guyer <guyer@nist.gov>
 # Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 # Author: James Warren   <jwarren@nist.gov>
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
 # History
 # 
 # modified   by  rev reason
 # ---------- --- --- -----------
 # 2006-10-31 JEG 1.0 original
 # 
 # ########################################################################
 ##

__docformat__ = 'restructuredtext'

class Stepper:
    def __init__(self, iterates=()):
        self.iterates = iterates
        
    def sweepFn(iterates, dt):
        residual = 0
        for var, eqn, bcs in iterates:
            residual = max(residual, eqn.sweep(var=var, dt=dt, boundaryConditions=bcs))
             
        return residual
    sweepFn = staticmethod(sweepFn)
         
    def successFn(iterates, dt, dtPrev, elapsed, *args, **kwargs):
        pass
    successFn = staticmethod(successFn)
         
    def failFn(iterates, dt, *args, **kwargs):
        pass
    failFn = staticmethod(failFn)

    def _lowerBound(self, dt):
        dt = max(dt, self.dtMin)
        if self.elapsed + dt == self.elapsed:
            raise "step size underflow: %g + %g == %g" % (self.elapsed, dt, self.elapsed)
            
        return dt
        
    def _step(self, dt, dtPrev, sweepFn, failFn, *args, **kwargs):
        sweepFn(iterates=self.iterates, dt=dt, *args, **kwargs) 
        return dt, dt
         
    def step(self, dt, dtTry=None, dtMin=None, dtPrev=None,
             sweepFn=None, successFn=None, failFn=None, *args, **kwargs):
        sweepFn = sweepFn or self.sweepFn
        successFn = successFn or self.successFn
        failFn = failFn or self.failFn
     
        dtTry = dtTry or dtMin or dt
        dtPrev = dtPrev or dtMin
        self.dtMin = dtMin or 0.
        
        self.elapsed = 0.
         
        while self.elapsed < dt:
            dtMax = dt - self.elapsed 
            if dtTry > dtMax:
                dtSave = dtTry
                dtTry = dtMax
            else:
                dtSave = None
            
            for var, eqn, bcs in self.iterates:
                var.updateOld()
                 
            dtPrev, dtTry = self._step(dt=dtTry, dtPrev=dtPrev,  
                                       sweepFn=sweepFn, failFn=failFn,
                                       *args, **kwargs)
                                      
            print "dtPrev:", dtPrev, "dtTry:", dtTry
            
            self.elapsed += dtPrev
                                
            successFn(iterates=self.iterates, 
                      dtPrev=dtPrev, elapsed=self.elapsed, dt=dt, *args, **kwargs)
                      
            dtTry = max(dtTry, self.dtMin)

        return dtSave or dtPrev, dtSave or dtTry
