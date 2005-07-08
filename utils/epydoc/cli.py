#!/usr/bin/env python
#
# objdoc: epydoc command-line interface
# Edward Loper
#
# Created [03/15/02 10:31 PM]
# $Id$
#

# Note: if you change this docstring, check that you didn't break
# _usage.
# Note: As it is, the usage message fits in an 80x24 window, but if it
# gets any bigger, it won't.
"""
Command-line interface for epydoc.

Usage::

 epydoc [OPTIONS] MODULES...
 
     MODULES...                The Python modules to document.
     --html                    Generate HTML output (default).
     --latex                   Generate LaTeX output.
     --pdf                     Generate pdf output, via LaTeX.
     --check                   Run documentation completeness checks.
     -o DIR, --output DIR      The output directory.
     -n NAME, --name NAME      The documented project's name.
     -u URL, --url URL         The documented project's url.
     -t PAGE, --top PAGE       The top page for the HTML documentation.
     -c SHEET, --css SHEET     CSS stylesheet for HTML files.
     --private-css SHEET       CSS stylesheet for private objects.
     --inheritance STYLE       The format for showing inherited objects.
     -V, --version             Print the version of epydoc.
     -h, -?, --help, --usage   Display this usage message.
     -h TOPIC, --help TOPIC    Display information about TOPIC (docformat,
                               css, inheritance, usage, or version).

 See the epydoc(1) man page for a complete list of options.

@var PROFILE: Whether or not to run the profiler.
@var TESTS: The lists of tests that can be run with
    C{'epydoc --check'}.

@var _encountered_internal_error: A global variable recording whether
any internal errors have been detected.  If this variable is set to
true, then L{cli} will issue a warning once it completes running.
"""
__docformat__ = 'epytext en'

##################################################
## Constants
##################################################

# Use "%(tex)s" to include the latex filename, "%(ps)s" for the
# postscript filename, etc.  (You must include the "s" after the close
# parenthasis).
LATEX_COMMAND = r"echo x | latex '\batchmode\input %(tex)s'"
MAKEINDEX_COMMAND = 'makeindex -q %(idx)s'
DVIPS_COMMAND = 'dvips -q %(dvi)s -o %(ps)s -G0 -Ppdf'
PS2PDF_COMMAND = ('ps2pdf -sPAPERSIZE=letter -dMaxSubsetPct=100 '+
                  '-dSubsetFonts=true -dCompatibilityLevel=1.2 '+
                  '-dEmbedAllFonts=true %(ps)s %(pdf)s')

# Testing (pdftex):
#LATEX_COMMAND = r"echo x | pdftex '\batchmode\input %(tex)s'"
#MAKEINDEX_COMMAND = 'makeindex -q %(idx)s'
#DVIPS_COMMAND = 'true'
#PS2PDF_COMMAND = 'true'

## This is a more verbose version of LATEX_COMMAND.
DEBUG_LATEX_COMMAND = r"echo x | latex %(tex)s"

PROFILE=0

# What tests can we run?
TESTS=('basic', 'types', 'vars', 'private', 'authors', 'versions', 'all')

##################################################
## Command-Line Interface
##################################################
import sys, os.path, re, getopt

# Include support for Zope, if it's available.
try: import ZODB
except: pass

def cli():
    """
    Command line interface for epydoc.
    
    @rtype: C{None}
    """
    # Parse the command line arguments.
    options = _parse_args()

    # Import all the specified modules.
    modules = _import(options['modules'], options['verbosity'])

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
    docmap = _make_docmap(modules, options)

    # Perform the requested action.
    if options['action'] == 'html': _html(docmap, options)
    elif options['action'] == 'check': _check(docmap, options)
    elif options['action'] == 'latex': _latex(docmap, options, 'latex')
    elif options['action'] == 'dvi': _latex(docmap, options, 'dvi')
    elif options['action'] == 'ps': _latex(docmap, options, 'ps')
    elif options['action'] == 'pdf': _latex(docmap, options, 'pdf')
    else: raise ValueError('Unknown action %r' % options['action'])

    # Report any internal errors.
    if _encountered_internal_error:
        estr = ("!! An internal error occured.  To see the exception "+
                "that caused the !!\n!! error, use the '--debug' "+
                "option.                                 !!")
        print >>sys.stderr, '\n'+'!'*70
        print >>sys.stderr, estr
        print >>sys.stderr, '!'*70+'\n'

_encountered_internal_error = 0
def _internal_error(e=None):
    """
    Print a warning message about an internal error.
    @return: The return value from calling C{func}
    """
    if isinstance(e, KeyboardInterrupt): raise
    global _encountered_internal_error
    _encountered_internal_error = 1
    if sys.stderr.softspace: print >>sys.stderr
    if e: print >>sys.stderr, "INTERNAL ERROR: %s" % e
    else: print >>sys.stderr, "INTERNAL ERROR"
        
def _usage(exit_code=1):
    """
    Display a usage message.

    @param exit_code: An exit status that will be passed to
        C{sys.exit}.
    @type exit_code: C{int}
    @rtype: C{None}
    """
    if exit_code == 0: stream = sys.stdout
    else: stream = sys.stderr
    NAME = os.path.split(sys.argv[0])[-1]
    if NAME == '(imported)': NAME = 'epydoc'
    usage = __doc__.split('Usage::\n')[-1].replace('epydoc', NAME)
    usage = re.sub(r'\n\s*@[\s\S]*', '', usage)
    usage = re.sub(r'\n ', '\n', usage)
    print >>stream, '\nUsage:', usage.strip()+'\n'
    sys.exit(exit_code)

def _usage_error(estr, exit_code=1):
    """
    Issue an error message, and exit.
    """
    progname = os.path.basename(sys.argv[0])
    if '\n' in estr:
        estr = '\n%s\nRun "%s -h" for usage.\n' % (estr.strip(), progname)
    else:
        estr = '%s; run "%s -h" for usage.' % (estr.strip(), progname)
        from epydoc.markup import wordwrap
        estr = '\n'+wordwrap(estr)
    print >>sys.stderr, estr
    sys.exit(exit_code)

def _help(arg):
    """
    Display a speficied help message, and exit.

    @param arg: The name of the help message to display.  Currently,
        only C{"css"} and C{"usage"} are recognized.
    @type arg: C{string}
    @rtype: C{None}
    """
    arg = arg.strip().lower()
    if arg == 'css':
        from epydoc.css import STYLESHEETS
        print '\nThe following built-in CSS stylesheets are available:'
        names = STYLESHEETS.keys()
        names.sort()
        maxlen = max(*[len(name) for name in names])
        format = '    %'+`-maxlen-1`+'s %s'
        for name in names:
            print format % (name, STYLESHEETS[name][1])
        print
    elif arg == 'version':
        _version()
    elif arg in ('checks', 'tests', 'test', 'check'):
        print '\nBy default, epydoc checks to make sure that all public'
        print 'objects have descriptions.  The following additional tests'
        print 'can be specified with the "--tests" option:'
        print '    - private: Check private objects.'
        print '    - vars: Check variables and parameters.'
        print '    - types: Check that variables and parameters have types.'
        print '    - authors: Check that all modules have authors.'
        print '    - versions: Check that all modules have versions.'
        print '    - all: Run all tests\n'
    elif arg in ('inheritance', 'inheritence'):
        print '\nThe following inheritance formats are currently supported:'
        print '    - grouped: inherited objects are gathered into groups,'
        print '      based on what class they were inherited from.'
        print '    - listed: inherited objects are listed in a short list'
        print '      at the end of their section.'
        print '    - included: inherited objects are mixed in with '
        print '      non-inherited objects.\n'
    elif arg in ('docformat', 'doc_format', 'doc-format'):
        print '\n__docformat__ is a module variable that specifies the markup'
        print 'language for the docstrings in a module.  Its value is a '
        print 'string, consisting the name of a markup language, optionally '
        print 'followed by a language code (such as "en" for English).  Epydoc'
        print 'currently recognizes the following markup language names:'
        import epydoc.objdoc
        for format in epydoc.objdoc.KNOWN_DOCFORMATS:
            print '  - %s' % format
        print
    else:
        _usage(0)
    sys.exit(0)
    
def _version():
    """
    Display the version information, and exit.
    
    @rtype: C{None}
    """
    import epydoc
    print "Epydoc version %s" % epydoc.__version__
    sys.exit(0)

def _check_css(cssname):
    """
    If C{cssname} is not valid, then issue an error and exit.
    """
    if cssname is None: return
    if os.path.isfile(cssname): return
    from epydoc.css import STYLESHEETS
    if STYLESHEETS.has_key(cssname): return

    # We couldn't find it.
    print >>sys.stderr, '\nError: CSS file %s not found\n' % cssname
    sys.exit(1)

def _parse_args():
    """
    Process the command line arguments; return a dictionary containing
    the relevant info.

    @return: A dictionary mapping from configuration parameters to
        values.  If a parameter is specified on the command line, then
        that value is used; otherwise, a default value is used.
        Currently, the following configuration parameters are set:
        C{target}, C{modules}, C{verbosity}, C{prj_name}, C{output},
        C{show_imports}, C{frames}, C{private}, C{debug}, C{top},
        C{list_classes_separately}, C{docformat}, C{inheritance},
        C{autogen_vars}, and C{test}.
    @rtype: C{None}
    """
    # Default values.
    options = {'target':None, 'modules':[], 'verbosity':1,
               'prj_name':'', 'action':'html', 'tests':{'basic':1},
               'show_imports':0, 'frames':1, 'private':None,
               'list_classes_separately': 0, 'debug':0,
               'docformat':None, 'top':None, 'inheritance': None,
               'ignore_param_mismatch': 0, 'alphabetical': 1}

    # Get the command-line arguments, using getopts.
    shortopts = 'c:h:n:o:t:u:Vvq?:'
    longopts = ('html latex dvi ps pdf check '+
                'output= name= url= top= css= private-css= private_css= '+
                'docformat= doc-format= doc_format= private '+
                'show_imports no_private no-private '+
                'builtins no-frames no_frames noframes debug '+
                'help= usage= helpfile= help-file= help_file= '+
                'separate-classes separate_classes '+
                'quiet show-imports show_imports '+
                'target= version verbose '+
                'navlink= nav_link= nav-link= '+
                'command-line-order command_line_order '+
                'inheritance= inheritence= '+
                'ignore_param_mismatch ignore-param-mismatch '+
                'test= tests= checks=').split()
    try:
        (opts, modules) = getopt.getopt(sys.argv[1:], shortopts, longopts)
    except getopt.GetoptError, e:
        if e.opt in ('h', '?', 'help', 'usage'): _usage(0)
        print >>sys.stderr, ('%s; run "%s -h" for usage' %
                              (e,os.path.basename(sys.argv[0])))
        sys.exit(1)

    # Parse the arguments.
    for (opt, val) in opts:
        if opt in ('--builtins',):
            modules += sys.builtin_module_names
            modules.remove('__main__')
        elif opt in ('--check',): options['action'] = 'check'
        elif opt in ('--command-line-order', '--command_line_order'):
            options['alphabetical'] = 0
        elif opt in ('--css', '-c'): options['css'] = val
        elif opt in ('--debug',):
            options['debug'] = 1
            options['verbosity'] += 4
        elif opt in ('--dvi',): options['action'] = 'dvi'
        elif opt in ('--docformat', '--doc-format', '--doc_format'):
            from epydoc.objdoc import set_default_docformat
            set_default_docformat(val)
        elif opt in ('--help', '-?', '--usage', '-h'): _help(val)
        elif opt in ('--helpfile', '--help-file', '--help_file'):
            options['help'] = val
        elif opt in ('--html',): options['action'] = 'html'
        elif opt in ('--ignore_param_mismatch', '--ignore-param-mismatch'):
            options['ignore_param_mismatch'] = 1
        elif opt in ('--inheritance', '--inheritence'):
            options['inheritance']=val.lower()
        elif opt in ('--latex',): options['action']='latex'
        elif opt in ('--name', '-n'): options['prj_name']=val
        elif opt in ('--navlink', '--nav-link', '--nav_link'):
            options['prj_link'] = val
        elif opt in ('--no-frames', '--no_frames', '--noframes'):
            options['frames'] = 0
        elif opt in ('--no-private', '--no_private'): options['private']=0
        elif opt in ('--output', '--target', '-o'): options['target']=val
        elif opt in ('--pdf',): options['action'] = 'pdf'
        elif opt in ('--private',):
            options['private'] = 1
            options['tests']['private'] = 1
        elif opt in ('--private-css', '--private_css'):
            options['private_css'] = val
        elif opt in ('--ps',): options['action'] = 'ps'
        elif opt in ('--quiet', '-q'): options['verbosity'] -= 1
        elif opt in ('--separate-classes', '--separate_classes'):
            options['list_classes_separately'] = 1
        elif opt in ('--show-imports', '--show_imports'):
            options['show_imports'] = 1
        elif opt in ('--test', '--tests', '--checks'):
            for test in re.split('\s*[,\s]\s*', val.lower()):
                options['tests'][test] = 1
        elif opt in ('-t', '--top'): options['top'] = val
        elif opt in ('--url', '-u'): options['prj_url']=val
        elif opt in ('--verbose', '-v'): options['verbosity'] += 1
        elif opt in ('--version', '-V'): _version()
        else:
            _usage()

    #//////////////////////////////
    # Default Values
    #//////////////////////////////
    # This section deals with default values that depend on what
    # action we're performing.

    # Pick a default target directory, if none was specified.
    if options['target'] is None:
        options['target'] = options['action']

    # Pick a default for private/no-private, if none was specified.
    if options['private'] is None:
        if options['action'] == 'html': options['private'] = 1
        else: options['private'] = 0

    # Pick a default for inheritance, if none was specified.
    if options['inheritance'] is None:
        if options['action'] == 'html': options['inheritance'] = 'grouped'
        else: options['inheritance'] = 'listed'

    #//////////////////////////////
    # Validity Checks
    #//////////////////////////////

    # Make sure inheritance has a valid value
    if options['inheritance'] not in ('grouped', 'listed', 'included'):
        estr = 'Bad inheritance style.  Valid options are '
        estr += 'grouped, listed, included'
        _usage_error(estr)

    # Make sure tests has a valid value.
    for test in options['tests'].keys():
        if test not in TESTS:
            estr = 'Bad epydoc test %r.  Valid tests are:' % test
            for t in TESTS: estr += '\n  - %s' % t
            _usage_error(estr)

    # Check that the options all preceed the filenames.
    for m in modules:
        if m == '-': break
        elif m[0:1] == '-':
            _usage_error('options must preceed modules')

    # Check the CSS file(s)
    _check_css(options.get('css'))
    _check_css(options.get('private_css'))
        
    # Make sure we got some modules.
    modules = [m for m in modules if m != '-']
    if len(modules) == 0:
        _usage_error('no modules specified')
    options['modules'] = modules

    return options

def _import(module_names, verbosity):
    """
    @return: A list of the modules contained in the given files.
        Duplicates are removed.  Order is preserved.
    @rtype: C{list} of C{module}
    @param module_names: The list of module filenames.
    @type module_names: C{list} of C{string}
    @param verbosity: Verbosity level for tracing output.
    @type verbosity: C{int}
    """
    from epydoc.imports import import_module, find_modules

    # First, expand packages.
    for name in module_names[:]:
        if os.path.isdir(name):
            # In-place replacement.
            index = module_names.index(name)
            new_modules = find_modules(name)
            if new_modules:
                module_names[index:index+1] = new_modules
            elif verbosity >= 0:
                if sys.stderr.softspace: print >>sys.stderr
                print  >>sys.stderr, 'Error: %r is not a pacakge' % name

    if verbosity > 0:
        print >>sys.stderr, 'Importing %s modules.' % len(module_names)
    modules = []
    progress = _Progress('Importing', verbosity, len(module_names))
    
    for name in module_names:
        progress.report(name)
        # Import the module, and add it to the list.
        try:
            module = import_module(name)
            if module not in modules: modules.append(module)
            elif verbosity > 2:
                if sys.stderr.softspace: print >>sys.stderr
                print >>sys.stderr, '  (duplicate)'
        except ImportError, e:
            if verbosity >= 0:
                if sys.stderr.softspace: print >>sys.stderr
                print  >>sys.stderr, e

    if len(modules) == 0:
        print >>sys.stderr, '\nError: no modules successfully loaded!'
        sys.exit(1)
    return modules

def _make_docmap(modules, options):
    """
    Construct the documentation map for the given modules.

    @param modules: The modules that should be documented.
    @type modules: C{list} of C{Module}
    @param options: Options from the command-line arguments.
    @type options: C{dict}
    """
    from epydoc.objdoc import DocMap, report_param_mismatches

    verbosity = options['verbosity']
    document_bases = 1
    document_autogen_vars = 1
    inheritance_groups = (options['inheritance'] == 'grouped')
    inherit_groups = (options['inheritance'] != 'grouped')
    d = DocMap(verbosity, document_bases, document_autogen_vars,
               inheritance_groups, inherit_groups)
    if options['verbosity'] > 0:
        print  >>sys.stderr, ('Building API documentation for %d modules.'
                              % len(modules))
    progress = _Progress('Building docs for', verbosity, len(modules))
    
    for module in modules:
        progress.report(module.__name__)
        # Add the module.  Catch any exceptions that get generated.
        try: d.add(module)
        except Exception, e:
            if options['debug']: raise
            else: _internal_error(e)
        except:   
            if options['debug']: raise
            else: _internal_error()

    if not options['ignore_param_mismatch']:
        if not report_param_mismatches(d):
            estr = '    (To supress these warnings, '
            estr += 'use --ignore-param-mismatch)'
            print >>sys.stderr, estr

    return d

def _run(cmd, options):
    from epydoc.markup import wordwrap
    if '|' in cmd: name = cmd.split('|')[1].strip().split(' ', 1)[0]
    else: name = cmd.strip().split(' ', 1)[0]
    if options['verbosity'] == 1:
        print >>sys.stderr, 'Running %s...' % name
    elif options['verbosity'] > 1:
        cmd_str = wordwrap(`cmd`, 10+len(name)).lstrip()
        print >>sys.stderr, 'Running %s' % cmd_str.rstrip()

    exitcode = os.system(cmd)
    if exitcode != 0:
        raise OSError('%s failed: exitcode=%s' % (name, exitcode))

def _latex(docmap, options, format):
    """
    Create the LaTeX documentation for the objects in the given
    documentation map.  

    @param docmap: A documentation map containing the documentation
        for the objects whose API documentation should be created.
    @param options: Options from the command-line arguments.
    @type options: C{dict}
    @param format: One of C{'latex'}, C{'dvi'}, C{'ps'}, or C{'pdf'}.
    """
    from epydoc.latex import LatexFormatter

    # Create the documenter, and figure out how many files it will
    # generate.
    latex_doc = LatexFormatter(docmap, **options)
    num_files = latex_doc.num_files()
        
    # Write documentation.
    if options['verbosity'] > 0:
        print  >>sys.stderr, ('Writing LaTeX docs (%d files) to %r.' %
                              (num_files, options['target']))
    progress = _Progress('Writing', options['verbosity'], num_files)
    latex_doc.write(options['target'], progress.report)

    # Run latex, makeindex, dvi, ps, and pdf, as appropriate.
    oldpath = os.path.abspath(os.curdir)
    try:
        try:
            # Filenames (used by the external commands)
            filenames = {'tex': 'api.tex', 'idx': 'api.idx',
                         'dvi': 'api.dvi', 'ps': 'api.ps',
                         'pdf': 'api.pdf'}
        
            # latex -> dvi
            if format in ('dvi', 'ps', 'pdf'):
                # Go into the output directory.
                os.chdir(options['target'])

                if options['debug']: latex_command = DEBUG_LATEX_COMMAND
                else: latex_command = LATEX_COMMAND
                
                _run(latex_command % filenames, options)
                _run(MAKEINDEX_COMMAND % filenames, options)
                _run(latex_command % filenames, options)
                _run(latex_command % filenames, options)
                
            # dvi -> postscript
            if format in ('ps', 'pdf'):
                _run(DVIPS_COMMAND % filenames, options)

            # postscript -> pdf
            if format in ('pdf',): 
                _run(PS2PDF_COMMAND % filenames, options)

        except OSError, e:
            print  >>sys.stderr, 'Error: %s' % e
            if not options['debug']:
                print 'Running epydoc with the --debug option may',
                print 'give more informative output.'
            sys.exit(1)
        except Exception, e:
            if options['debug']: raise
            else: _internal_error(e)
        except:   
            if options['debug']: raise
            else: _internal_error()
    finally:
        os.chdir(oldpath)

def _html(docmap, options):
    """
    Create the HTML documentation for the objects in the given
    documentation map.  

    @param docmap: A documentation map containing the documentation
        for the objects whose API documentation should be created.
    @param options: Options from the command-line arguments.
    @type options: C{dict}
    """
    from epydoc.html import HTMLFormatter

    # Create the documenter, and figure out how many files it will
    # generate.
    html_doc = HTMLFormatter(docmap, **options)
    num_files = html_doc.num_files()

    # Write documentation.
    if options['verbosity'] > 0:
        print  >>sys.stderr, ('Writing HTML docs (%d files) to %r.' %
                              (num_files, options['target']))
    progress = _Progress('Writing', options['verbosity'], num_files, 1)
    try: html_doc.write(options['target'], progress.report)
    except OSError, e:
        print >>sys.stderr, '\nError writing docs:\n%s\n' % e
    except IOError, e:
        print >>sys.stderr, '\nError writing docs:\n%s\n' % e
    except Exception, e:
        if options['debug']: raise
        else: _internal_error(e)
    except:   
        if options['debug']: raise
        else: _internal_error()

def _check(docmap, options):
    """
    Run completeness checks on the objects in the given documentation
    map.  By default, C{_check} checks for docstrings in all public
    modules, classes, functions, and properties.  Additional checks
    can be added with the C{'tests'} option:

      - C{private}: Also checks private objects.
      - C{vars}: Also checks variables, parameters, and return values.

    @param docmap: A documentation map containing the documentation
        for the objects whose API documentation should be created.
    @param options: Options from the command-line arguments.
    @type options: C{dict}
    """
    from epydoc.checker import DocChecker
    
    # Run completeness checks.
    if options['verbosity'] > 0:
        print  >>sys.stderr, 'Performing completeness checks...'
    checker = DocChecker(docmap)

    if options['tests'].get('all'):
        for test in TESTS: options['tests'][test] = 1

    # Run the checks
    checks = 0
    if (options['tests'].get('basic') or
        options['tests'].get('vars') or
        options['tests'].get('private')):
        checks |= (DocChecker.MODULE | DocChecker.CLASS |
                   DocChecker.FUNC | DocChecker.PROPERTY |
                   DocChecker.DESCR_LAZY | DocChecker.PUBLIC)
    if options['tests'].get('private'): checks |= DocChecker.PRIVATE
    if options['tests'].get('vars'): checks |= DocChecker.ALL_T
    if options['tests'].get('types'):
        checks |= DocChecker.ALL_T
        DocChecker.TYPE
    passed_checks = checker.check(checks)

    if options['tests'].get('authors'):
        checks = DocChecker.MODULE | DocChecker.PUBLIC | DocChecker.AUTHOR
        if options['tests'].get('private'): checks |= DocChecker.PRIVATE
        passed_checks = checker.check(checks) and passed_checks
    
    if options['tests'].get('versions'):
        checks = DocChecker.MODULE | DocChecker.PUBLIC | DocChecker.VERSION
        if options['tests'].get('private'): checks |= DocChecker.PRIVATE
        passed_checks = checker.check(checks) and passed_checks
    
    if passed_checks and options['verbosity'] > 0:
        print >>sys.stderr, '  All checks passed!'
           
class _Progress:
    """

    The progress meter that is used by C{cli} to report its progress.
    It prints the status to C{stderrr}.  Depending on the verbosity,
    setting it will produce different outputs.

    To update the progress meter, call C{report} with the name of the
    object that is about to be processed.
    """
    def __init__(self, action, verbosity, total_items, html_file=0):
        """
        Create a new progress meter.

        @param action: A string indicating what action is performed on
            each objcet.  Examples are C{"writing"} and C{"building
            docs for"}.
        @param verbosity: The verbosity level.  This controls what the
            progress meter output looks like.
        @param total_items: The total number of items that will be
            processed with this progress meter.  This is used to let
            the user know how much progress epydoc has made.
        @param html_file: Whether to assume that arguments are html
            file names, and munge them appropriately.
        """
        self._action = action
        self._verbosity = verbosity
        self._total_items = total_items
        self._item_num = 1
        self._html_file = 0

    def report(self, argument):
        """
        Update the progress meter.
        @param argument: The object that is about to be processed.
        """
        if self._verbosity <= 0: return
        
        if self._verbosity==1:
            if self._item_num == 1 and self._total_items <= 70:
                sys.stderr.write('  [')
            if (self._item_num % 60) == 1 and self._total_items > 70:
                sys.stderr.write('  [%3d%%] ' %
                                 (100.0*self._item_num/self._total_items))
            sys.stderr.write('.')
            sys.stderr.softspace = 1
            if (self._item_num % 60) == 0 and self._total_items > 70:
                print >>sys.stderr
            if self._item_num == self._total_items:
                if self._total_items <= 70: sys.stderr.write(']')
                print >>sys.stderr
        elif self._verbosity>1:
            TRACE_FORMAT = (('  [%%%dd/%d]' % (len(`self._total_items`),
                                               self._total_items))+
                            ' %s %%s' % self._action)

            if self._html_file:
                (dir, file) = os.path.split(argument)
                (root, d) = os.path.split(dir)
                if d in ('public', 'private'):
                    argument = os.path.join(d, file)
                else:
                    fname = argument
            
            print >>sys.stderr, TRACE_FORMAT % (self._item_num, argument)
        self._item_num += 1
        
if __name__ == '__main__':
    if PROFILE:
        import profile
        profile.run('cli()', '/tmp/profile.out')
        import pstats
        p = pstats.Stats('/tmp/profile.out')
        p.strip_dirs().sort_stats('time', 'cum').print_stats(60)
        p.strip_dirs().sort_stats('cum', 'time').print_stats(60)
    else:
        cli()
