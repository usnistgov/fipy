#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "dump.py"
 #                                    created: 1/10/04 {10:23:17 AM} 
 #                                last update: 1/16/04 {12:16:56 PM} 
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

import cPickle
import os

def write(data, fileName):
    fileStream = os.popen('gzip --fast -c > %s.gz' % (fileName), 'w')
    cPickle.dump(data, fileStream, 0)
    fileStream.close()

def read(fileName):
    fileStream = os.popen('gunzip --fast -c < %s.gz' % (fileName), 'r')
    data = cPickle.load(fileStream)
    fileStream.close()
    return data

    
