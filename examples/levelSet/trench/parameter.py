#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "parameter.py"
 #                                    created: 8/19/04 {10:29:10 AM} 
 #                                last update: 8/19/04 {4:01:07 PM} { 1:23:41 PM}
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

"""

This file keeps all the geometry and material properties necessary to
run a trench level set electrochemistry problem in one
location. Everything is in S.I. units.

"""

materialProperties = {
    "Faraday's constant" : 9.6e4,
    "gas constant" : 8.314
    }

geometryParameters = {
    "depth of trench" : 0.5e-6,
    "trench aspect ratio" : 2,
    "trench spacing" : 0.4e-6,
    "boundary layer depth" : 0.3e-6
    }

metalParameters = {
    "atomic volume" : 7.1e-6,  
    "change" : 2,
    "bulk molar concentration" : 278,
    "diffusion coefficient" : 5.6e-10,
    }

experimentalParameters = {
    "exchange current density" : {
        "constant" : 0.26,
        "accelerator dependence" : 45},

    "transfer coefficient" : {
        "constant" : 0.5,
        "accelerator dependence" : 0},
    
    "temperature" : 298,
    "overpotential" : -0.3,
    "initial accelerator coverage" : 0.1
    }

controlParameters = {
    "duration of simulation" : 100,
    "cfl number" : 0.5,
    "cell size" : 0.5e-7,
    "number of cells in narrow band" : 10,
    "cells below trench" : 10
    }
