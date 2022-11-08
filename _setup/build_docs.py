from distutils.core import Command
import os

__all__ = ["build_docs"]

class build_docs(Command):

    description = "build the FiPy documentation"

    # List of option tuples: long name, short name (None if no short
    # name), and help string.
    user_options = [('pdf', None, "compile the PDF variant of the documentation"),
                    ('html', None, "compile the HTML variant of the documentation"),
                    ('cathartic', None, "rewrite all the files (default is to only rewrite changed files)"),
                   ]

    def initialize_options (self):
        self.pdf = 0
        self.html = 0
        self.cathartic = 0

    def finalize_options (self):
        pass

    def run (self):
        import sphinx.cmd.build
        import sphinx.ext.apidoc
        
        sphinx_args = ['-P', '-n', '-c', 'documentation/', '.']
        apidoc_args = []
        
        if self.cathartic:
            sphinx_args = ['-a', '-E'] + sphinx_args
            apidoc_args = ['--force'] + apidoc_args
            
        sphinx.ext.apidoc.main(['--output-dir=fipy/generated', '--suffix=rst']
                    + apidoc_args + ['fipy'])
        sphinx.ext.apidoc.main(['--output-dir=documentation/tutorial/package/generated', '--suffix=rst']
                    + apidoc_args + ['documentation/tutorial/package'])

        if self.html:
            sphinx.cmd.build.main(['-b', 'redirecting_html'] + sphinx_args + ['documentation/_build/html/'])

        if self.pdf:
            try:
                sphinx.cmd.build.main(['-b', 'latex'] + sphinx_args + ['documentation/_build/latex/'])
            except SystemExit:
                pass
            
            outdir = os.path.join('documentation', '_build', 'latex')
            
            from docutils.core import publish_file

            for xtra in ("LICENSE", "DISCLAIMER"):
                publish_file(source_path="%s.rst" % xtra,
                             destination_path=os.path.join(outdir, "%s.tex" % xtra),
                             reader_name='standalone',
                             parser_name='restructuredtext',
                             writer_name='latex',
                             settings_overrides= {
                                 'template': 'documentation/_templates/empty.tex'
                             })

            savedir = os.getcwd()
            
            os.chdir(outdir)
                
            os.system("pdflatex fipy")
            os.system("pdflatex fipy")
            os.system("pdflatex fipy")
            os.system("makeindex -s python.ist fipy")
            os.system("makeindex -s python.ist modfipy")
            os.system("pdflatex fipy")
            os.system("pdflatex fipy")
                
            os.chdir(savedir)
