# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build
RST2TEXTXT    = $(SOURCEDIR)/LICENSE.rst $(SOURCEDIR)/DISCLAIMER.rst
TEMPLATE      = "$(SOURCEDIR)/_templates/empty.tex"
R2LOPTS       = --use-latex-citations --new-column-widths

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Generate isolated TeX from ReST for \input in nist.sty
%.tex.txt: %.rst
	rst2latex.py --template=$(TEMPLATE) $(R2LOPTS) $< > $@

# Break circular dependency
%.rst: ;

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile $(RST2TEXTXT:.rst=.tex.txt)
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
