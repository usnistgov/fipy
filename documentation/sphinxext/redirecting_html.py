from docutils import nodes
import urlparse

from sphinx.builders.html import StandaloneHTMLBuilder

class RedirectingHTMLBuilder(StandaloneHTMLBuilder):
    name = 'redirecting_html'

    def write_doc(self, docname, doctree):
        """Replace URLs external to NIST with a redirection

        Based slightly on `CheckExternalLinksBuilder.check()`
        """
        for node in doctree.traverse(nodes.reference):
            try:
                uri = node['refuri']
                uri = urlparse.urlparse(uri)
                if uri.scheme in ["http", "https"]:
                    if not uri.netloc.endswith("nist.gov"):
                        node['refuri'] = "/cgi-bin/redirect.py?url=" + uri.geturl()
            except KeyError:
                continue

        super(RedirectingHTMLBuilder, self).write_doc(docname, doctree)

def setup(app):
    app.add_builder(RedirectingHTMLBuilder)
