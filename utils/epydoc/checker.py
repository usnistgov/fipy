#
# objdoc: epydoc documentation completeness checker
# Edward Loper
#
# Created [01/30/01 05:18 PM]
# $Id$
#

"""
Documentation completeness checker.  This module defines a single
class, C{DocChecker}, which can be used to check the that specified
classes of objects are documented.
"""
__docformat__ = 'epytext en'

##################################################
## Imports
##################################################

import re, sys, os.path, string
from xml.dom.minidom import Text as _Text
from types import ModuleType as _ModuleType
from types import FunctionType as _FunctionType
from types import BuiltinFunctionType as _BuiltinFunctionType
from types import BuiltinMethodType as _BuiltinMethodType
from types import MethodType as _MethodType
from types import StringType as _StringType

from epydoc.uid import UID, Link
from epydoc.objdoc import ModuleDoc, FuncDoc, PropertyDoc
from epydoc.objdoc import ClassDoc, Var, Raise, ObjDoc

# The following methods may be undocumented:
_NO_DOCS = ['__hash__', '__repr__', '__str__', '__cmp__']

# The following methods never need descriptions, authors, or
# versions:
_NO_BASIC = ['__hash__', '__repr__', '__str__', '__cmp__']

# The following methods never need return values:
_NO_RETURN = ['__init__', '__hash__', '__repr__', '__str__',
              '__cmp__']

# The following methods don't need parameters documented:
_NO_PARAM = ['__cmp__']

def _is_private(str):
    """
    @return: True if C{str} is the name of a public object.
    @rtype: C{boolean}
    @param str: The name to check.
    @type str: C{string}
    """
    if str == '...': return 0
    for piece in str.split('.'):
        if piece[0] == '_' and piece[-1] != '_':
            return 1
    return 0

class DocChecker:
    """
    Documentation completeness checker.  C{DocChecker} can be used to
    check that specified classes of objects are documented.  To check
    the documentation for a group of objects, you should create a
    C{DocChecker} from a L{DocMap<objdoc.DocMap>} that documents those
    objects; and then use the L{check} method to run specified checks
    on the objects' documentation.

    What checks are run, and what objects they are run on, are
    specified by the constants defined by C{DocChecker}.  These
    constants are divided into three groups.  

      - Type specifiers indicate what type of objects should be
        checked: L{MODULE}; L{CLASS}; L{FUNC}; L{VAR}; L{IVAR};
        L{CVAR}; L{PARAM}; and L{RETURN}.
      - Public/private specifiers indicate whether public or private
        objects should be checked: L{PUBLIC} and L{PRIVATE}.
      - Check specifiers indicate what checks should be run on the
        objects: L{TYPE}; L{DESCR}; L{DESCR_LAZY}; L{AUTHOR};
        and L{VERSION}.

    The L{check} method is used to perform a check on the
    documentation.  Its parameter is formed by or-ing together at
    least one value from each specifier group:

        >>> checker.check(DocChecker.MODULE | DocChecker.PUBLIC |
        ...               DocChecker.DESCR)
        
    To specify multiple values from a single group, simply or their
    values together:
    
        >>> checker.check(DocChecker.MODULE | DocChecker.CLASS |
        ...               DocChecker.FUNC | DocChecker.DESCR_LAZY |
        ...               DocChecker.PUBLIC)

    @group Types: MODULE, CLASS, FUNC, VAR, IVAR, CVAR, PARAM,
        RETURN, ALL_T
    @type MODULE: C{int}
    @cvar MODULE: Type specifier that indicates that the documentation
        of modules should be checked.
    @type CLASS: C{int}
    @cvar CLASS: Type specifier that indicates that the documentation
        of classes should be checked.
    @type FUNC: C{int}
    @cvar FUNC: Type specifier that indicates that the documentation
        of functions should be checked.
    @type VAR: C{int}
    @cvar VAR: Type specifier that indicates that the documentation
        of module variables should be checked.
    @type IVAR: C{int}
    @cvar IVAR: Type specifier that indicates that the documentation
        of instance variables should be checked.
    @type CVAR: C{int}
    @cvar CVAR: Type specifier that indicates that the documentation
        of class variables should be checked.
    @type PARAM: C{int}
    @cvar PARAM: Type specifier that indicates that the documentation
        of function and method parameters should be checked.
    @type RETURN: C{int}
    @cvar RETURN: Type specifier that indicates that the documentation
        of return values should be checked.
    @type ALL_T: C{int}
    @cvar ALL_T: Type specifier that indicates that the documentation
        of all objects should be checked.

    @group Checks: TYPE, AUTHOR, VERSION, DESCR_LAZY, DESCR, ALL_C
    @type TYPE: C{int}
    @cvar TYPE: Check specifier that indicates that every variable and
        parameter should have a C{@type} field.
    @type AUTHOR: C{int}
    @cvar AUTHOR: Check specifier that indicates that every object
        should have an C{author} field.
    @type VERSION: C{int}
    @cvar VERSION: Check specifier that indicates that every object
        should have a C{version} field.
    @type DESCR_LAZY: C{int}
    @cvar DESCR_LAZY: Check specifier that indicates that every object
        should have a description.  However, it is permissible for
        functions and methods that have a C{@return} field to not have
        a description, since a description may be generated from the
        C{@return} field.
    @type DESCR: C{int}
    @cvar DESCR: Check specifier that indicates that every object
        should have a description.  
    @type ALL_C: C{int}
    @cvar ALL_C: Check specifier that indicates that  all checks
        should be run.

    @group Publicity: PUBLIC, PRIVATE, ALL_P
    @type PUBLIC: C{int}
    @cvar PUBLIC: Specifier that indicates that public objects should
        be checked.
    @type PRIVATE: C{int}
    @cvar PRIVATE: Specifier that indicates that private objects should
        be checked.
    @type ALL_P: C{int}
    @cvar ALL_P: Specifier that indicates that both public and private
        objects should be checked.
    """
    # Types
    MODULE = 1
    CLASS  = 2
    FUNC   = 4
    VAR    = 8
    IVAR   = 16
    CVAR   = 32
    PARAM  = 64
    RETURN = 128
    PROPERTY = 256
    ALL_T  = 1+2+4+8+16+32+64+128+256

    # Checks
    TYPE = 256
    AUTHOR = 1024
    VERSION = 2048
    DESCR_LAZY = 4096
    _DESCR_STRICT = 8192
    DESCR = DESCR_LAZY + _DESCR_STRICT
    ALL_C = 256+512+1024+2048+4096+8192

    # Private/public
    PRIVATE = 16384
    PUBLIC = 32768
    ALL_P = PRIVATE + PUBLIC

    ALL = ALL_T + ALL_C + ALL_P

    def __init__(self, docmap):
        """
        Create a new C{DocChecker} that can be used to run checks on
        the documentation of the objects documented by C{docmap}

        @param docmap: A documentation map containing the
            documentation for the objects to be checked.
        @type docmap: L{DocMap<objdoc.DocMap>}
        """
        self._docmap = docmap
        
        # Sort & filter the documentation from the docmap.
        docs = docmap.items()
        docs.sort()
        self._docs = [d for (u,d) in docs if u.is_module() or
                docmap.has_key(u.module())]

        # Initialize instance variables
        self._checks = 0
        self._last_warn = None
        self._out = sys.stdout
        self._num_warnings = 0

    def check(self, checks = None):
        """
        Run the specified checks on the documentation of the objects
        contained by this C{DocChecker}'s C{DocMap}.  Any errors found
        are printed to standard out.

        @param checks: The checks that should be run on the
            documentation.  This value is constructed by or-ing
            together the specifiers that indicate which objects should
            be checked, and which checks should be run.  See the
            L{module description<checker>} for more information.
            If no checks are specified, then a default set of checks
            will be run.
        @type checks: C{int}
        @return: True if no problems were found.
        @rtype: C{boolean}
        """
        if checks == None:
            return (self.check(DocChecker.MODULE | DocChecker.CLASS |
                               DocChecker.FUNC | DocChecker.DESCR_LAZY |
                               DocChecker.PUBLIC) and
                    self.check(DocChecker.PARAM | DocChecker.VAR |
                               DocChecker.IVAR | DocChecker.CVAR |
                               DocChecker.RETURN | DocChecker.DESCR |
                               DocChecker.TYPE | DocChecker.PUBLIC))

        self._checks = checks
        self._last_warn = None
        for doc in self._docs:
            uid = doc.uid()
            if isinstance(doc, ModuleDoc):
                self._check_module(doc)
            elif isinstance(doc, ClassDoc):
                self._check_class(doc)
            elif isinstance(doc, FuncDoc):
                self._check_func(doc)
            elif isinstance(doc, PropertyDoc):
                self._check_property(doc)
            else:
                raise AssertionError("Don't know how to check%r" % doc)
        if self._last_warn is not None: print
        return (self._last_warn is None)

    def _warn(self, error, name):
        """
        Print a warning about the object named C{name}.

        @param error: The error message to print.
        @type error: C{string}
        @param name: The name of the object that generated the
            warning.
        @type name: C{string}
        @rtype: C{None}
        """
        name = str(name)
        if name != self._last_warn:
            self._num_warnings += 1
            if self._last_warn is not None: self._out.write('\n')
            self._out.write(name + '.'*max(1,60-len(name)) + error)
            self._last_warn = name
        else:
            self._out.write(', %s' % error)

    def _check_name_publicity(self, name):
        """
        @return: True if an object named C{name} should be checked,
            given the public/private specifiers.
        @rtype: C{boolean}
        @param name: The name of the object to check.
        @type name: C{string}
        """
        if (_is_private(name) and
            not (self._checks & DocChecker.PRIVATE)): return 0
        if (not _is_private(name) and
            not (self._checks & DocChecker.PUBLIC)): return 0
        return 1

    def _check_basic(self, doc):
        """
        Check the description, author, version, and see-also fields of
        C{doc}.  This is used as a helper function by L{_check_module},
        L{_check_class}, and L{_check_func}.

        @param doc: The documentation that should be checked.
        @type doc: L{ObjDoc}
        @rtype: C{None}
        """
        if (self._checks & DocChecker.DESCR) and (not doc.descr()):
            if ((self._checks & DocChecker._DESCR_STRICT) or
                (not isinstance(doc, FuncDoc)) or
                (not doc.returns().descr())):
                self._warn('No descr', doc.uid().name())
        if self._checks & DocChecker.AUTHOR:
            for field in doc.fields():
                if 'author' in field.tags: break
            else:
                self._warn('No authors', doc.uid().name())
        if self._checks & DocChecker.VERSION:
            for field in doc.fields():
                if 'version' in field.tags: break
            else:
                self._warn('No version', doc.uid().name())
            
    def _check_module(self, doc):
        """
        Run checks on the module whose UID is C{doc}.
        
        @param doc: The UID of the module to check.
        @type doc: L{UID}
        @rtype: C{None}
        """
        if not self._check_name_publicity(doc.uid().name()): return
        if self._checks & DocChecker.MODULE:
            self._check_basic(doc)
        if self._checks & DocChecker.VAR:
            for v in doc.variables():
                self._check_var(v, doc.uid().name())
        
    def _check_class(self, doc):
        """
        Run checks on the class whose UID is C{doc}.
        
        @param doc: The UID of the class to check.
        @type doc: L{UID}
        @rtype: C{None}
        """
        if not self._check_name_publicity(doc.uid().name()): return
        if self._checks & DocChecker.CLASS:
            self._check_basic(doc)
        if self._checks & DocChecker.IVAR:
            for v in doc.ivariables():
                self._check_var(v, doc.uid().name())
        if self._checks & DocChecker.CVAR:
            for v in doc.cvariables():
                self._check_var(v, doc.uid().name())

    def _check_property(self, doc):
        if not self._check_name_publicity(doc.uid().name()): return
        if self._checks & DocChecker.PROPERTY:
            self._check_basic(doc)

    def _check_var(self, var, name, check_type=1):
        """
        Run checks on the variable whose documentation is C{var} and
        whose name is C{name}.
        
        @param var: The documentation for the variable to check.
        @type var: L{Var}
        @param name: The name of the variable to check.
        @type name: C{string}
        @param check_type: Whether or not the variable's type should
            be checked.  This is used to allow varargs and keyword
            parameters to have no type specified.
        @rtype: C{None}
        """
        if not self._check_name_publicity(name): return
        if var == None: return
        if not self._check_name_publicity(var.name()): return
        if var.name() == 'return':
            if (var.type() and
                var.type().to_plaintext(None).strip().lower() == 'none'):
                return
        if (self._checks & DocChecker.DESCR) and (not var.descr()):
            self._warn('No descr', name+'.'+var.name())
        if ((self._checks & DocChecker.TYPE) and (not var.type()) and
            check_type):
            self._warn('No type', name+'.'+var.name())
            
    def _check_func(self, doc):
        """
        Run checks on the function whose UID is C{doc}.
        
        @param doc: The UID of the function to check.
        @type doc: L{UID}
        @rtype: C{None}
        """
        name = doc.uid().name()
        if not self._check_name_publicity(name): return
        if (doc.uid().is_any_method() and 
            doc != self._docmap.documented_ancestor(doc.uid())): return
        if (self._checks & DocChecker.FUNC and
            not doc.has_docstring() and
            doc.uid().shortname() not in _NO_DOCS):
            self._warn('No docs', name)
            return
        if (self._checks & DocChecker.FUNC and
            doc.uid().shortname() not in _NO_BASIC):
                self._check_basic(doc)
        if (self._checks & DocChecker.RETURN and
            doc.uid().shortname() not in _NO_RETURN):
                self._check_var(doc.returns(), name)
        if (self._checks & DocChecker.PARAM and
            doc.uid().shortname() not in _NO_PARAM):
            if doc.uid().is_method():
                for v in doc.parameters()[1:]:
                    self._check_var(v, name)
            else:
                for v in doc.parameters():
                    self._check_var(v, name)
            self._check_var(doc.vararg(), name, 0)
            self._check_var(doc.kwarg(), name, 0)
