## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "driver.py"
 #                                    created: 3/17/05 {11:55:06 AM} 
 #                                last update: 3/18/05 {10:18:14 AM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  epydocker.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 #  Description: 
 #    
 #    Replacement for the 'epydoc' command-line tool.
 #       Driven from a Python setup.py script (or other), rather than from a shell
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2005-03-17 JEG 1.0 original
 # ###################################################################
 ##

import sys

import epydoc.cli

_options = {'target':None, 'modules':[], 'verbosity':1,
           'prj_name':'', 'action':'html', 'tests':{'basic':1},
           'show_imports':0, 'frames':1, 'private': 0,
           'list_classes_separately': 0, 'debug':0,
           'docformat':None, 'top':None, 'inheritance': 'listed',
           'ignore_param_mismatch': 0, 'alphabetical': 1}

def epylatex(module_names, options = {}):
    for key in _options.keys():
        if not options.has_key(key):
            options[key] = _options[key]
    
    modules = epydoc.cli._import(module_names, options['verbosity'])
               
    # Record the order of the modules in options.
    from epydoc.uid import make_uid
    options['modules'] = muids = []
    for m in modules:
        try:
            muids.append(make_uid(m))
        except:
            if sys.stderr.softspace: print >>sys.stderr
            print >>sys.stderr, 'Failed to create a UID for %s' % m

    # Build their documentation
    docmap = epydoc.cli._make_docmap(modules, options)
    
    _latex(docmap, options)
    
    if epydoc.cli._encountered_internal_error:
        estr = ("!! An internal error occured.  To see the exception "+
                "that caused the !!\n!! error, use the '--debug' "+
                "option.                                 !!")
        print >>sys.stderr, '\n'+'!'*70
        print >>sys.stderr, estr
        print >>sys.stderr, '!'*70+'\n'

def _latex(docmap, options):
    """
    Create the LaTeX documentation for the objects in the given
    documentation map.  

    @param docmap: A documentation map containing the documentation
        for the objects whose API documentation should be created.
    @param options: Options from the command-line arguments.
    @type options: C{dict}
    """
    from latex import LatexFormatter

    # Create the documenter, and figure out how many files it will
    # generate.
    latex_doc = LatexFormatter(docmap, **options)
    num_files = latex_doc.num_files()
        
    # Write documentation.
    if options['verbosity'] > 0:
        print  >>sys.stderr, ('Writing LaTeX docs (%d files) to %r.' %
                              (num_files, options['target']))
    progress = epydoc.cli._Progress('Writing', options['verbosity'], num_files)
    latex_doc.write(options['target'], progress.report)
