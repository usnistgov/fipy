#!/usr/bin/env python
#
# import: module import support for epydoc
# Edward Loper
#
# Created [10/06/02 02:14 AM]
# $Id$
#

"""
Module import support for epydoc.
"""
__docformat__ = 'epytext en'

import sys, re, types, os.path

_VALID_MODULE_NAME = re.compile(r'^[a-zA-Z_]\w*(\.[a-zA-Z_]\w*)*$')

def find_modules(dirname):
    """
    If C{dirname} contains a package, then return the filenames for
    the package, its modules, and all of its descendant packages and
    modules.  If C{dirname} does not contain a package, then return an
    empty list.

    @rtype: C{list} of C{string}
    @return: The filenames for the package in C{dirname}, its modules,
        and all of its descendant packages and modules.
    @type dirname: C{string}
    @param dirname: The directory name to search.
    """
    if not os.path.isdir(dirname): return []

    # Find and import the __init__.py file.  We need to do this to get
    # its __path__ variable.
    try:
        if os.path.isfile(os.path.join(dirname, '__init__.py')):
            pkg = import_module(os.path.join(dirname, '__init__.py'))
        elif os.path.isfile(os.path.join(dirname, '__init__.pyw')):
            pkg = import_module(os.path.join(dirname, '__init__.pyw'))
        elif os.path.isfile(os.path.join(dirname, '__init__.pyc')):
            pkg = import_module(os.path.join(dirname, '__init__.pyc'))
        else:
            return [] # Not a package
    except:
        return [] # Error importing package

    modules = _find_modules(pkg).keys()
    modules.sort()
    return [pkg.__name__] + modules

def _find_modules(pkg):
    """
    Helper function for L{find_modules}.
    """
    pkgname = pkg.__name__
    
    # Get a list of all files/directories on the package's path.
    paths = []
    for dir in pkg.__path__:
        if os.path.isdir(dir):
            paths += [os.path.join(dir, file) for
                      file in os.listdir(dir)]

    # Search the files/directories for modules & subdirs.
    modules = {}  # a set (dict from name to [1,0])
    subdirectories = []
    for filepath in paths:
        file = os.path.split(filepath)[-1]
        if os.path.isdir(filepath):
            if re.match(r'\w+', file):
                subpkgname = '%s.%s' % (pkg.__name__, file)
                try: subpkg = _import_module(subpkgname)
                except: continue
                if not hasattr(subpkg, '__path__'): continue
                modules[subpkgname] = 1
                modules.update(_find_modules(subpkg))
        elif not re.match(r'\w+.py.?', file):
            continue # Ignore things like ".#foo.py" or "a-b.py"
        elif file == '__init__.py' or file[:-1] == '__init__.py':
            continue
        elif file[-3:] == '.py':
            modules['%s.%s' % (pkgname, file[:-3])] = 1
        elif file[-4:-1] == '.py':
            modules['%s.%s' % (pkgname, file[:-4])] = 1
        
    return modules

def import_module(name_or_filename):
    """
    @return: The module with the given module name or filename.
    
        - If C{name_or_filename} is the name of a module (such as
          C{os.path}), then return that module.
        - If C{name_or_filename} is the name of a file (such as
          C{epytext.py} or C{multiarray.so}), then return the module
          defined by that file.

    If C{name_or_filename} names both a file and a module, then
    C{import_module} will treat it as a file name.
    
    @rtype: C{list} of C{module}
    @type name_or_filename: C{string}
    @param name_or_filename: The module name or filename specifying
        which module to import.
    @raise ImportError: If there was any problem importing the given
        object.  In particular, an C{ImportError} is raised if the
        given file does not exist; if the given file name does not
        name a valid module; or if importing the module causes any
        exception to be raised.
    """
    # If we've already imported it, then just return it.
    if sys.modules.has_key(name_or_filename):
        return sys.modules[name_or_filename]
    
    # Try importing it as a file name.
    if os.path.exists(name_or_filename):
        old_sys_path = sys.path[:]
        (basedir, name) = _find_module_from_filename(name_or_filename)
        sys.path.insert(0, basedir)
        try: return _import_module(name)
        finally: sys.path = old_sys_path

    # Try importing it as a module name.
    if _VALID_MODULE_NAME.match(name_or_filename):
        return _import_module(name_or_filename)

    # Report an error.
    raise ImportError('Error importing %r: ' % name_or_filename +
                      'bad module name or file not found.')

def _import_module(name):
    """
    @return: the module with the given name.  C{import_module} makes
        some attempts to prevent the imported module from modifying
        C{sys} and C{__builtins__}.  In the future, more sandboxing
        might be added (e.g., using the C{rexec} module).
    @rtype: C{module}
    @param name: The name of the module to import.  C{name} is a
        fully qualified module name, such as C{os.path}.
    @type name: C{string}
    @raise ImportError: If there was any problem importing the given
        object.  In particular, if importing the module causes any
        exception to be raised.
    """
    # If we've already imported it, then just return it.
    if sys.modules.has_key(name): return sys.modules[name]

    # Save some important values.  This helps prevent the module that
    # we import from breaking things *too* badly.  Note that we do
    # *not* save sys.modules.
    old_sys = sys.__dict__.copy()
    old_sys_path = sys.path[:]
    if type(__builtins__) == types.DictionaryType:
        old_builtins = __builtins__.copy()
    else:
        old_builtins = __builtins__.__dict__.copy()

    # Supress input and output.  (These get restored when we restore
    # sys to old_sys).  
    sys.stdin = sys.stdout = sys.stderr = _dev_null
    sys.__stdin__ = sys.__stdout__ = sys.__stderr__ = _dev_null

    # Remove any command-line arguments
    sys.argv = ['(imported)']

    try:
        # Make sure that we have a valid name
        if not _VALID_MODULE_NAME.match(name):
            raise ImportError('Error importing %r: bad module name' % name)
            
        # Import the module.  Note that if "name" has a package
        # component, then this just gives us the top-level object.
        try:
            topLevel = __import__(name)
        except KeyboardInterrupt:
            raise # don't capture keyboard interrupts!
        except Exception, e:
            raise ImportError('Error importing %r: %s' % (name, e))
        except SystemExit, e:
            raise ImportError('Error importing %r: %s' % (name, e))
        except:
            raise ImportError('Error importing %r')
        
        # If "name" has a package component, then we have to manually
        # go down the package tree.
        pieces = name.split(".")[1:]
        m = topLevel
        for p in pieces:
            try: m = getattr(m, p)
            except AttributeError:
                estr = 'Error importing %r: getattr failed' % name
                raise ImportError(estr)
        return m
    finally:
        # Restore the important values that we saved.
        if type(__builtins__) == types.DictionaryType:
            __builtins__.clear()
            __builtins__.update(old_builtins)
        else:
            __builtins__.__dict__.clear()
            __builtins__.__dict__.update(old_builtins)
        sys.__dict__.clear()
        sys.__dict__.update(old_sys)
        sys.path = old_sys_path

class _DevNull:
    """
    A "file-like" object that discards anything that is written and
    always reports end-of-file when read.  C{_DevNull} is used by
    L{import_module} to discard output when importing modules; and to
    ensure that stdin appears closed.
    """
    def __init__(self):
        self.closed = 1
        self.mode = 'r+'
        self.softspace = 0
        self.name='</dev/null>'
    def close(self): pass
    def flush(self): pass
    def read(self, size=0): return ''
    def readline(self, size=0): return ''
    def readlines(self, sizehint=0): return []
    def seek(self, offset, whence=0): pass
    def tell(self): return 0L
    def truncate(self, size=0): pass
    def write(self, str): pass
    def writelines(self, sequence): pass
    xreadlines = readlines
_dev_null = _DevNull()

def _find_module_from_filename(filename):
    """
    Break a module/package filename into a base directory and a module
    name.  C{_find_module_from_filename} checks directories in the
    filename to see if they contain C{"__init__.py"} files; if they
    do, then it assumes that the module is part of a package, and
    returns the full module name.  For example, if C{filename} is
    C{"/tmp/epydoc/imports.py"}, and the file
    C{"/tmp/epydoc/__init__.py"} exists, then the base directory will
    be C{"/tmp/"} and the module name will be C{"epydoc.imports"}.
    
    @return: A pair C{(basedir, module)}, where C{basedir} is the base
        directory from which the module can be imported; and C{module}
        is the name of the module itself.    
    @rtype: C{(string, string)}

    @param filename: The filename that contains the module.
        C{filename} can be a directory (for a package); a C{.py} file;
        a C{.pyc} file; a C{.pyo} file; or an C{.so} file.
    @type filename: C{string}
    """
    # Normalize the filename
    filename = os.path.normpath(os.path.abspath(filename))

    # Split the file into (basedir, module, ext), and check the extension.
    (basedir, file) = os.path.split(filename)
    (module, ext) = os.path.splitext(file)
    if not (ext[-3:] == '.py' or ext[-4:-1] == '.py' or
            ext[-3:] == '.so'):
        raise ImportError('Error importing %r: ' % filename +
                          'not a Python module')

    # Is it a package?
    if module == '__init__':
        (basedir, module) = os.path.split(basedir)
    
    # If there's a package, then find its base directory.
    if (os.path.exists(os.path.join(basedir, '__init__.py')) or
        os.path.exists(os.path.join(basedir, '__init__.pyc')) or
        os.path.exists(os.path.join(basedir, '__init__.pyw'))):
        package = []
        while os.path.exists(os.path.join(basedir, '__init__.py')):
            (basedir,dir) = os.path.split(basedir)
            if dir == '': break
            package.append(dir)
        package.reverse()
        module = '.'.join(package+[module])

    return (basedir, module)

if __name__ == '__main__':
    # A few quick tests.
    print import_module('os.path')
    print import_module('distutils')
    print import_module('/usr/lib/python2.1/linecache.py')
    print import_module('/usr/lib/python2.1/distutils/cmd.py')
    print import_module('/usr/lib/python2.1/distutils/__init__.py')
    #print import_module('/usr/lib/python2.1/distutils/')
    print import_module('/usr/lib/python2.1/site-packages/'+
                        'Numeric/multiarray.so')

    
