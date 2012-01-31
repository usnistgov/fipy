""" sphinx extension to implement bibtex blocks

Implements two new directives::

    .. biblisted:: myrefs.bib
        :style: default

        ref1
        ref2
        ref3

which results in the generation of 3 citation targets with the references
extracted from myrefs.bib and formatted according to the `default` style.

.. bibmissing:: myrefs.bib
    :style: default

This directive takes no content, and watches for missing citations during the
build.  If it has a citation for the missing reference, it inserts the citation
during build.

The option ``sort`` to either directive::

    .. bibmissing:: myrefs.bib
        :sort:

results in the citations found being sorted according to the specified or
default ``style``.  See the ``bibstuff`` code and documentation for how to
specify styles.

Again, for both directives, you can specify more than one bib file separated by
commas::

    .. bibmissing:: myrefs.bib, morerefs.bib

All heavy lifting done by a slighly modified version of the ``bibstuff``
library, by Alan Isaac.

:copyright: Copyright 2010 by Matthew Brett
:license: simplified BSD
"""
DEBUG = False
if DEBUG:
    # This is ipython 0.10 (not 0.11, for which the interface has changed)
    from IPython.Shell import IPShellEmbed
    ipshell = IPShellEmbed([])

import os
from os import path
import sys
import codecs

from docutils.parsers.rst import nodes, directives, Parser
from docutils.utils import new_document
from docutils.nodes import fully_normalize_name

from sphinx import addnodes
from sphinx.util.compat import Directive
from sphinx.util.nodes import make_refnode

from .. import bibfile, bibgrammar, bibstyles

class CiteMaker(object):
    """ Class to make citations from bibliographies

    It has to be possible to pickle this class because it will be saved with the
    doctree.

    Parameters
    ----------
    bib_str : str
        string with bibtex bibliography text
    style_str : str, optional
        name of the citation style
    extra_styles : None or mapping, optional
        name, module URI mapping for any extra styles.
    """
    def __init__(self, bib_str, style_str='default', extra_styles=None):
        self.bib_str = bib_str
        self.style_str = style_str
        if extra_styles is None:
            extra_styles = {}
        self.extra_styles = extra_styles
        self._bibobj = None
        self._citation_manager = None

    @property
    def bibobj(self):
        if self._bibobj is None:
            # create object to store parsed .bib file
            parsed_bibfile = bibfile.BibFile()
            # store parsed .bib files in the bibfile_processor
            bibgrammar.Parse(self.bib_str, parsed_bibfile)
            self._bibobj = parsed_bibfile
        return self._bibobj

    @property
    def citation_manager(self):
        if not self._citation_manager is None:
            return self._citation_manager
        # Get style; we do this to avoid storing the module in the class, so it
        # can be pickled
        style_module = self.style_getter(self.style_str)
        # Override any fancy reference semantics for link target
        class CM(style_module.CitationManager):
            def get_citation_label(self,entry,citation_template=None):
                return '.. [' + entry.citekey + ']\n'
        self._citation_manager = CM([self.bibobj])
        return self._citation_manager

    def style_getter(self, name):
        # Fill cs dictionary on the fly to allow class to be pickled
        cs = custom_styles(self.extra_styles)
        if name in cs:
            return cs[name]
        return bibstyles.from_name(name)

    def entry_for_ref(self, ref):
        return self.bibobj.get_entry_by_citekey(ref)

    def citestr_for_ref(self, ref):
        entry = self.entry_for_ref(ref)
        if entry is None:
            return None
        return self.citestr_for_entry(entry)

    def citestr_for_entry(self, entry):
        cm = self.citation_manager
        return cm.format_citation(entry)


def custom_styles(styles_named):
    """ Read style modules from name, module string mapping
    """
    styles = {}
    for name, mod in styles_named.items():
        if not hasattr(mod, 'CitationManager'):
            __import__(mod)
            mod = sys.modules[mod]
        styles[name] = mod
    return styles


# Copied from hg sphinx environment.py
def relfn2path(self, filename, docname=None):
    """Return paths to a file referenced from a document, relative to
    documentation root and absolute.

    Absolute filenames are relative to the source dir, while relative
    filenames are relative to the dir of the containing document.
    """
    if filename.startswith('/') or filename.startswith(os.sep):
        rel_fn = filename[1:]
    else:
        docdir = path.dirname(self.doc2path(docname or self.docname,
                                            base=None))
        rel_fn = path.join(docdir, filename)
    try:
        return rel_fn, path.join(self.srcdir, rel_fn)
    except UnicodeDecodeError:
        # the source directory is a bytestring with non-ASCII characters;
        # let's try to encode the rel_fn in the file system encoding
        enc_rel_fn = rel_fn.encode(sys.getfilesystemencoding())
        return rel_fn, path.join(self.srcdir, enc_rel_fn)


def arg_to_filenames(arg):
    """ Parse input argument line to contents

    Filenames are separated by commas
    """
    return [name.strip() for name in arg.split(',')]


def resolve_citations(doc, citations):
    """
    Link citations to/from their references.
    """
    for citation in citations:
        for label in citation['names']:
            if label in doc.citation_refs:
                reflist = doc.citation_refs[label]
                resolve_references(citation, reflist)


def resolve_references(note, reflist):
# # #     assert len(note['ids']) == 1
    id = note['ids'][0]
    for ref in reflist:
        if ref.resolved:
            continue
        ref.delattr('refname')
        ref['refid'] = id
        assert len(ref['ids']) == 1
        note.add_backref(ref['ids'][0])
        ref.resolved = 1
    note.resolved = 1


def rst2nodes(text, settings):
    new_doc = new_document('temp-string', settings)
    parser = Parser()
    parser.parse(text, new_doc)
    return new_doc.children


def make_citation(label, text, settings):
    name = fully_normalize_name(label)
    citation = nodes.citation(text)
    citation += nodes.label('', label)
    new_doc = new_document('temp-string', settings)
    parser = Parser()
    parser.parse(text, new_doc)
    citation['names'].append(name)
    citation += new_doc.children
    return citation


class BibListedDirective(Directive):

    has_content = True
    # Take arguments as one long string, we will parse into comma-separated
    # filenames later
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    option_spec = {
        'style': directives.unchanged,
        'sort': directives.flag,
        'encoding': directives.encoding,
    }

    def run(self):
        """ Extract references from directive content and insert citations

        Here our job is very easy, because we just replace the directive content
        with the citations that it has promised, and no more need be seen of
        this directive after that is done
        """
        document = self.state.document
        cite_maker, sort_flag, message = self.get_bibs_opts()
        if cite_maker is None:
            return [document.reporter.warning(message, line=self.lineno)]
        # Refs are just not-blank lines that don't start with #
        refs = []
        for line in self.content:
            ref = line.strip()
            if ref != '' and not ref.startswith('#'):
                refs.append(ref)
        # Find corresponding refs
        entries = []
        warnings = []
        for ref in refs:
            entry = cite_maker.entry_for_ref(ref)
            if entry is None:
                warnings.append(document.reporter.warning(
                    'cannot find ref %s' % ref,
                    line=self.lineno))
                continue
            entries.append(entry)
        cm = cite_maker.citation_manager
        if sort_flag:
            entries.sort(key=cm.sortkey)
        cite_strs = [cm.format_citation(entry) for entry in entries]
        cite_text = '\n'.join(cite_strs)
        own_doc = self.state.document
        citations = rst2nodes(cite_text, own_doc.settings)
        # Make backrefs to references
        resolve_citations(own_doc, citations)
        # We replace ourselves by the citations we provide
        return citations + warnings

    def get_bibs_opts(self):
        """ Read passed bib files into object, return with options

        This part is common between directives for bibtex references

        It just reads the bibtex files with the given encoding, makes an object
        to contain them.

        Returns
        -------
        cite_maker : ``CiteMaker`` instance, or None
            object wrapping read bibliography, or None if there was an error
            creating the instance
        sort_flag : bool
            If True, ``sort`` option was set
        message : str
            If cite_maker is None, error message explaining why.
        """
        settings = self.state.document.settings
        fnames = arg_to_filenames(self.arguments[0])
        style_str = self.options.get('style', 'default')
        sort_flag = 'sort' in self.options
        encoding = self.options.get('encoding', settings.input_encoding)
        # Read bibfiles into string
        env = settings.env
        bib_str = ''
        for fname in fnames:
            rel_filename, filename = relfn2path(env, fname)
            env.note_dependency(rel_filename)
            if not os.path.exists(filename):
                return (None,
                        False,
                        'bibtex file %s does not exist' % filename)
            fp = codecs.open(filename, 'r', encoding)
            try:
                bibcontents = fp.read()
            except UnicodeDecodeError:
                return (None,
                        False,
                        'bibtex file %s should be %s encoding' %
                        (filename, encoding))
            finally:
                fp.close()
            bib_str += bibcontents
        extra_styles = env.config.bibref_styles
        cite_maker = CiteMaker(bib_str, style_str, extra_styles)
        return cite_maker, sort_flag, ''


class BibrefProvider(nodes.General, nodes.Element):
    """ Placeholder node for something that can provide references
    """
    pass


class BibMissingDirective(BibListedDirective):
    """ Directive to get missing references from bibtex files

    Now our job is much harder, compared to the case where the references are
    specified, because we have to defer providing the citations until after all
    the documents have been read, and we have a list of citations that don't
    otherwise exist in the doctree.
    """
    has_content = False

    def run(self):
        """ Initialize provider cache, return placeholder

        Placeholder holds the place for the eventual set of citations.
        """
        document = self.state.document
        cite_maker, sort_flag, message = self.get_bibs_opts()
        if cite_maker is None:
            return [document.reporter.warning(message, line=self.lineno)]
        env = document.settings.env
        # Id for this provider, pcache pair.  We need this id because the saving
        # and loading of the doctree means that a hash won't correspond
        prov_id = 'bibref-provider-%s' % env.new_serialno('bibref-provider')
        # Provider cache
        if not hasattr(env, 'bibref_providers'):
            env.bibref_providers = []
        pcache = dict(
            id = prov_id,
            cite_maker = cite_maker,
            docname = env.docname,
            document = document,
            sort_flag = sort_flag,
            citations = [],
            claimed = {})
        env.bibref_providers.append(pcache)
        provider = BibrefProvider('', ids=[prov_id])
        return [provider]


def env_purge_doc(app, env, docname):
    """ Remove cache of provider information in environment

    The providers live a dual life, with a placeholder ``BibrefProvider`` node,
    and the real information cached in the environment.  Here we've been asked
    to purge stuff for the document `docname`.  We remove any cached values for
    providers in this document.
    """
    if not hasattr(env, 'bibref_providers'):
        return
    env.bibref_providers = [provider for provider in env.bibref_providers
                            if provider['docname'] != docname]


def doctree_read(app, doctree):
    """ Give each pending citation its docname

    We need to do this, because later we have to make a link between this
    citation reference and citation that doesn't exist yet in the doctree.  To
    do this, we need the docnames for citation reference (here) and for the
    source of the (yet to be created) reference.
    """
    for pref in doctree.traverse(addnodes.pending_xref):
        if pref['reftype'] != 'citation':
            continue
        pref['docname'] = app.env.docname


def missing_reference(app, env, node, contnode):
    """ Check if any provider can provide this missing reference

    This event fires for every missing reference including missing citations.
    For missing citations, we ask each of the providers in turn whether they can
    provide this reference.  If so, we get a promised id for the citation,
    create a suitable reference node, and return it.  When we resolve the
    doctrees, we actually insert these promised citations (not here).  If we
    can't find a source for this reference, we return None, as the contract for
    ``missing-reference`` requires.
    """
    if not node['reftype'] == 'citation':
        return
    if not hasattr(env, 'bibref_providers'):
        return
    pcaches = env.bibref_providers
    if len(pcaches) == 0:
        return
    ref = node['reftarget']
    citation = None
    # Claimed?
    for pcache in pcaches:
        citation = pcache['claimed'].get(ref)
        if not citation is None:
            break
    if citation is None: # Not claimed
        # Search for it in available providers
        for pcache in env.bibref_providers:
            cite_maker = pcache['cite_maker']
            entry = cite_maker.entry_for_ref(ref)
            if entry is None:
                continue
            cite_txt = cite_maker.citestr_for_entry(entry)
            # Make a new citation node for this text
            own_doc = pcache['document']
            citation = rst2nodes(cite_txt, own_doc.settings)[0]
            citation['docname'] = pcache['docname']
            citation['entry'] = entry
            pcache['citations'].append(citation)
            # claim it
            pcache['claimed'][ref] = citation
            break
    if citation is None:
        # failed to find citation for this ref
        return # return None
    # Citation promised (and cached); make reference node for it
    fromdocname = node['docname'] # where ref is coming from
    docname = citation['docname'] # document it's pointing to
    labelid = citation['ids'][0]
    return make_refnode(app.builder, fromdocname, docname,
                        labelid, contnode)


def doctree_resolved(app, doctree, fromdocname):
    """ Put the promised references into provider trees

    During the ``missing-reference`` phase, the ``provider`` instances made some
    promises about the references they would provide.  Here we take those
    promises and fulfill them, by replacing the provider placeholders with the
    list of citations that they promised.

    If we asked for the references to be sorted, they'll be sorted according to
    the reference style, otherwise they come out in the order they are cited
    (whatever order the sphinx build process gave us).

    We try to make back references for the citations, but this can be difficult
    because we need to decide which of the possible references to back reference
    to.
    """
    env = app.env
    if not hasattr(env, 'bibref_providers'):
        return
    for bib_prov in doctree.traverse(BibrefProvider):
        id = bib_prov['ids'][0]
        pcache = [pc for pc in env.bibref_providers if pc['id'] == id][0]
        citations = pcache['citations']
        # resolve back references - just in this document for now
        resolve_citations(doctree.document, citations)
        # sort if requested
        if pcache['sort_flag']:
            cm = pcache['cite_maker'].citation_manager
            citations.sort(key=lambda x: cm.sortkey(x['entry']))
        bib_prov.replace_self(citations)


def setup(app):
    app.add_config_value('bibref_styles', {}, False)
    app.add_node(BibrefProvider)
    app.add_directive('biblisted', BibListedDirective)
    app.add_directive('bibmissing', BibMissingDirective)
    app.connect('env-purge-doc', env_purge_doc)
    app.connect('doctree-read', doctree_read)
    app.connect('missing-reference', missing_reference)
    app.connect('doctree-resolved', doctree_resolved)
