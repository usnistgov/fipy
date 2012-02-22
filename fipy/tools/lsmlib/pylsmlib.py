import pyximport
import numpy as np
pyximport.install(setup_args = {'options' :
                                {'build_ext' :
                                 {'libraries' : ['lsm_serial', 'lsm_toolbox'],
                                  'include_dirs' : np.get_include()
                                  }}})


from lsmlib import computeDistanceFunction2d_ as computeDistanceFunction2d
from lsmlib import computeExtensionFields2d_ as computeExtensionFields2d
from lsmlib import solveEikonalEquation2d_ as solveEikonalEquation2d
