#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "parser.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 9/3/04 {10:41:38 PM}
 # Stolen from:
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

"""

Parser routine to parse arguments for example files. Must use
'--blah=slomehting' otherwise will not be ablt to prune arguments'

"""

import optparse

import sys

def parse(larg, action = None, type = None, default = None):

    sarg = None
    tmpparser = optparse.OptionParser(option_list = [
        optparse.make_option(sarg, larg, action = action, type = type, dest = 'dest', default = default)],
                                      conflict_handler = 'ignore')
    
##    optparse.make_option('-e', '--numberOfElements', action = 'store', type = 'int', dest = 'Nele', default = numberOfElements),
##    optparse.make_option('-n', '--numberOfSteps', action = 'store', type = 'int', dest = 'steps', default = numberOfSteps),
##    optparse.make_option('-p', '--quiet', action = 'store_true', dest = 'quiet', default = False),
##    optparse.make_option('-i', '--inline', action = 'store_true', dest = 'inline', default = False)])

    sysargs = []
    for arg in sys.argv:
        if larg in arg:
            sysargs.append(arg)

    (options, args) = tmpparser.parse_args(sysargs)

    return options.dest

