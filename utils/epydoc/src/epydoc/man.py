#!/usr/bin/python2.2
#
# epydoc.py: manpage-style text output
# Edward Loper
#
# Created [01/30/01 05:18 PM]
# $Id$
#

"""
Documentation formatter that produces man-style documentation.

@note: This module is under development.  It generates incomplete
documentation pages, and is not yet incorperated into epydoc's
command-line interface.
"""
__docformat__ = 'epytext en'

##################################################
## Imports
##################################################

# system imports
import sys, xml.dom.minidom

# epydoc imports
import epydoc
from epydoc.uid import UID, Link, findUID, make_uid
from epydoc.imports import import_module
from epydoc.objdoc import DocMap, ModuleDoc, FuncDoc
from epydoc.objdoc import ClassDoc, Var, Raise, ObjDoc

##################################################
## Documentation -> Text Conversion
##################################################

class ManFormatter:
    def __init__(self, docmap, **kwargs):
        self._docmap = docmap

    #////////////////////////////////////////////////////////////
    # Basic Doc Pages
    #////////////////////////////////////////////////////////////
    
    def documentation(self, uid):
        if not self._docmap.has_key(uid):
            print '**NO DOCS ON %s **' % uid
            return
        doc = self._docmap[uid]
        
        if uid.is_module(): return self._modulepage(uid, doc)
        elif uid.is_class(): return self._classpage(uid, doc)
        elif uid.is_routine(): return self._routinepage(uid, doc)
        elif uid.is_variable(): return self._varpage(uid, doc)

    def _modulepage(self, uid, doc):
        str = self._name(uid)
        str += self._descr(uid, doc)
        str += self._funclist(doc.functions(), doc, 'FUNCTIONS')
        return str

    def _classpage(self, uid, doc):
        str = self._name(uid)
        str += self._descr(uid, doc)
        str += self._funclist(doc.methods(), doc, 'METHODS')
        str += self._funclist(doc.staticmethods(), doc, 'STATIC METHODS')
        str += self._funclist(doc.classmethods(), doc, 'CLASS METHODS')
        return str

    def _routinepage(self, uid, doc):
        str = self._name(uid)
        str += self._descr(uid, doc)
        return str

    def _varpage(self, uid, doc):
        str = self._name(uid)
        str += self._descr(uid, doc)
        return str

    #////////////////////////////////////////////////////////////
    # Functions
    #////////////////////////////////////////////////////////////
    
    def _funclist(self, functions, cls, title='FUNCTIONS'):
        str = self._title(title)
        numfuncs = 0
        
        for link in functions:
            fname = link.name()
            func = link.target()

            if func.is_method():
                container = func.cls()
                inherit = (container != cls.uid())
            else:
                inherit = 0
                try: container = func.module()
                except TypeError: container = None
            if not self._docmap.has_key(func):
                continue

            # If we don't have documentation for the function, then we
            # can't say anything about it.
            if not self._docmap.has_key(func): continue
            fdoc = self._docmap[func]

            # What does this method override?
            foverrides = fdoc.overrides()

            # Try to find a documented ancestor.
            inhdoc = self._docmap.documented_ancestor(func) or fdoc
            inherit_docs = (inhdoc is not fdoc)

            numfuncs += 1
            str += '    %s\n' % self._func_signature(self._bold(fname), fdoc)

            # Use the inherited docs for everything but the signature.
            fdoc = inhdoc

            fdescr=fdoc.descr()
            fparam = fdoc.parameter_list()[:]
            freturn = fdoc.returns()
            fraises = fdoc.raises()
            
            # Don't list parameters that don't have any extra info.
            f = lambda p:p.descr() or p.type()
            fparam = filter(f, fparam)

            # Description
            if fdescr:
                fdescr_str = fdescr.to_plaintext(None, indent=8)
                if fdescr_str.strip(): str += fdescr_str

            # Parameters
            if fparam:
                str += '        Parameters:\n'
                for param in fparam:
                    pname = param.name()
                    str += '            ' + pname
                    if param.descr():
                        pdescr = param.descr().to_plaintext(None, indent=12)
                        str += ' - %s' % pdescr.strip()
                    str += '\n'
                    if param.type():
                        ptype = param.type().to_plaintext(none, indent=16)
                        str += ' '*16+'(type=%s)\n' % ptype.strip()

            # Returns
            if freturn.descr():
                fdescr = freturn.descr().to_plaintext(None, indent=12)
                str += '        Returns:\n%s' % fdescr
                    
            if freturn.type():
                ftype = freturn.type().to_plaintext(None, indent=12)
                str += ("        Return Type: %s" % ftype.lstrip())

            ## Raises
            #if fraises:
            #    str += '        Raises:\n'
            #    for fraise in fraises:
            #        str += '      '
            #        str += ''+fraise.name()+' -\n'
            #        str += epytext.to_plaintext(fraise.descr(), 12)

            ## Overrides
            #if foverrides:
            #    str += '    <dl><dt><b>Overrides:</b></dt>\n'
            #    str += '      <dd>'+self._uid_to_href(foverrides)
            #    if inherit_docs:
            #        str += ' <i>(inherited documentation)</i>\n'
            #    str += '</dd>\n    </dl>\n'

        if numfuncs == 0: return ''

        return str
            
    def _func_signature(self, fname, fdoc, show_defaults=1):
        str = fname
        str += '('
        str += self._params_to_text(fdoc.parameters(), show_defaults)
        
        if fdoc.vararg():
            vararg_name = fdoc.vararg().name()
            if vararg_name != '...': vararg_name = '*%s' % vararg_name
            str += '%s, ' % vararg_name
        if fdoc.kwarg():
            str += '**%s, ' % fdoc.kwarg().name()
        if str[-1] != '(': str = str[:-2]

        return str + ')'
    
    def _params_to_text(self, parameters, show_defaults):
        str = ''
        for param in parameters:
            if type(param) in (type([]), type(())):
                sublist = self._params_to_text(param, 
                                               show_defaults)
                str += '(%s), ' % sublist[:-2]
            else:
                str += param.name()
                if show_defaults and param.default() is not None:
                    default = param.default()
                    if len(default) > 60:
                        default = default[:57]+'...'
                    str += '=%s' % default
                str += ', '
        return str

    #////////////////////////////////////////////////////////////
    # Helpers
    #////////////////////////////////////////////////////////////
    
    def _bold(self, text):
        """Format a string in bold by overstriking."""
        return ''.join([ch+'\b'+ch for ch in text])

    def _title(self, text):
        return '%s\n' % self._bold(text)

    def _kind(self, uid):
        if uid.is_package(): return 'package'
        elif uid.is_module(): return 'module'
        elif uid.is_class(): return 'class'
        elif uid.is_method() or uid.is_builtin_method(): return 'method'
        elif uid.is_routine(): return 'function'
        elif uid.is_variable(): return 'variable'
        else: raise AssertionError, 'Bad UID type for _name'

    def _name(self, uid):
        if uid.parent():
            parent = uid.parent()
            name = '%s %s in %s %s' % (self._kind(uid),
                                       self._bold(uid.shortname()),
                                       self._kind(parent),
                                       self._bold(parent.name()))
        else:
            name = '%s %s' % (self._kind(uid), self._bold(uid.name()))
        return '%s    %s\n\n' % (self._title('NAME'), name)

    def _descr(self, uid, doc):
        if not doc.descr(): return ''
        descr = doc.descr().to_plaintext(None, indent=4)
        return '%s%s' % (self._title('DESCRIPTION'), descr)

if __name__ == '__main__':
    docmap = DocMap(document_bases=1)

    uids = [findUID(name) for name in sys.argv[1:]]
    uids = [uid for uid in uids if uid is not None]
    for uid in uids: docmap.add(uid.value())

    formatter = ManFormatter(docmap)
    for uid in uids:
        print formatter.documentation(uid)
