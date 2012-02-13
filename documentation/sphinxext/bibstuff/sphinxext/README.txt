#######################
Sphinx bibtex extension
#######################

The extension is part of (my) Matthew Brett's edits to bibstuff.

The current code for my version is at http://github.com/matthew-brett/bibstuff .

To use the extension, make sure that this version of bibstuff is installed onto
your python path somewhere, and then add ``bibstuff.sphinxext.bibref`` to your
list of extensions in your sphinx ``conf.py`` file - for example, something like
this::

    extensions = ['sphinx.ext.pngmath', 'bibstuff.sphinxext.bibref']

There's an example of bibref in use at http://matthew.dynevor.org . The page
build code is here: http://gitorious.org/personal-www/dyneweb

Please see the ``bibref.py`` file docstring for more details.
