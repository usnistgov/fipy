############################################################
##  epydoc Makefile
##
##  Edward Loper
############################################################

##//////////////////////////////////////////////////////////////////////
## Configuration variables
##//////////////////////////////////////////////////////////////////////

# Python source files (don't include src/epydoc/test)
PY_SRC = src/epydoc/
EXAMPLES_SRC = $(wildcard doc/*.py)
DOCS = $(wildcard doc/*.html) $(wildcard doc/*.css) $(wildcard doc/*.png)

# What version of python to use?
PYTHON = python2.2

# The location of the webpage.
HOST = shell.sf.net
DIR = /home/groups/e/ep/epydoc/htdocs

# The current version of epydoc.
VERSION = $(shell python -c 'import epydoc; print epydoc.__version__')

# Base output directories
WEBDIR        = webpage
LATEX         = latex
HTML          = html

# Output subdirectories
HTML_API      = $(HTML)/api
HTML_EXAMPLES = $(HTML)/examples
HTML_STDLIB   = $(HTML)/stdlib
LATEX_API     = $(LATEX)/api
LATEX_STDLIB  = $(LATEX)/stdlib

EPYDOC = $(PYTHON) src/epydoc/cli.py
export PYTHONPATH=src/

##//////////////////////////////////////////////////////////////////////
## Usage
##//////////////////////////////////////////////////////////////////////

.PHONY: all usage clean distributions web webpage xfer local
.PHONY: checkdocs api-html api-pdf examples stdlib-html stdlib-pdf
.PHONY: test tests

all: usage
usage:
	@echo "Usage:"
	@echo "  make webpage -- build the webpage and copy it to sourceforge"
	@echo "    make api-html -- build the HTML docs for epydoc"
	@echo "    make api-pdf -- build the PDF docs for epydoc"
	@echo "    make examples -- build example API docs for the webpage"
	@echo "    make stdlib-html -- build HTML docs for the Python Standard Library"
	@echo "  make checkdocs -- check the documentation completeness"
	@echo "  make distributions -- build the distributions"
	@echo "  make clean -- remove all built files"
	@echo "  make test -- run unit tests"

##//////////////////////////////////////////////////////////////////////
## Clean
##//////////////////////////////////////////////////////////////////////

clean:
	$(MAKE) -C src clean
	rm -rf $(WEBDIR) $(HTML) $(LATEX)
	rm -rf .*.up2date

##//////////////////////////////////////////////////////////////////////
## Distributions
##//////////////////////////////////////////////////////////////////////

distributions: .distributions.up2date
.distributions.up2date: test $(PY_SRC) .webpage.up2date $(DOCS)
	$(MAKE) -C src distributions
	touch .distributions.up2date

##//////////////////////////////////////////////////////////////////////
## Web page
##//////////////////////////////////////////////////////////////////////

web: xfer
webpage: xfer
xfer: test .webpage.up2date stdlib-html
	rsync -arzv -e ssh $(WEBDIR)/* $(HOST):$(DIR)
	rsync -arzv -e ssh $(HTML_STDLIB)/ $(HOST):$(DIR)/stdlib

local: .webpage.up2date
	cp -r $(WEBDIR)/* /var/www/epydoc

checkdoc: checkdocs
checkdocs:
	epydoc --check --tests=vars,types $(PY_SRC)

.webpage.up2date: .api-html.up2date .examples.up2date .api-pdf.up2date \
		doc/epydoc-man.html doc/epydocgui-man.html $(DOCS)
	rm -rf $(WEBDIR)
	mkdir -p $(WEBDIR)
	cp -r $(DOCS) $(WEBDIR)
	cp -r $(HTML_API) $(WEBDIR)/api
	cp -r $(HTML_EXAMPLES) $(WEBDIR)/examples
	cp $(LATEX_API)/api.pdf $(WEBDIR)/epydoc.pdf
	touch .webpage.up2date

# Use plaintext docformat by default.  But this is overridden by the
# __docformat__ strings in each epydoc module.  (So just
# xml.dom.minidom and a few Docutils modules get plaintext
# docstrings).
api-html: .api-html.up2date
.api-html.up2date: $(PY_SRC)
	rm -rf $(HTML_API)
	mkdir -p $(HTML_API)
	$(EPYDOC) \
	       -o $(HTML_API) -n epydoc -u http://epydoc.sourceforge.net \
	       --inheritance=listed --navlink "epydoc 2.0"\
	       --css blue --private-css green -v --debug \
	       --docformat plaintext $(PY_SRC)
	touch .api-html.up2date

api-pdf: .api-pdf.up2date
.api-pdf.up2date: $(PY_SRC)
	rm -rf $(LATEX_API)
	mkdir -p $(LATEX_API)
	$(EPYDOC) --pdf -o $(LATEX_API) \
	       -n "Epydoc $(VERSION)" $(PY_SRC)
	touch .api-pdf.up2date

examples: .examples.up2date
.examples.up2date: $(EXAMPLES_SRC) $(PY_SRC)
	rm -rf $(HTML_EXAMPLES)
	mkdir -p $(HTML_EXAMPLES)
	$(EPYDOC) \
	       -o $(HTML_EXAMPLES) -n epydoc -u http://epydoc.sourceforge.net \
	       --no-private --css blue -t example --docformat=plaintext \
	       --navlink 'epydoc examples' doc/epydoc_example.py sre
	$(EPYDOC) -o $(HTML_EXAMPLES)/grouped \
	       --inheritance=grouped \
	       -n epydoc -u http://epydoc.sourceforge.net \
	       --no-private --css blue --debug \
	       --navlink 'epydoc examples' doc/inh_example.py
	$(EPYDOC) -o $(HTML_EXAMPLES)/listed \
	       --inheritance=listed \
	       -n epydoc -u http://epydoc.sourceforge.net \
	       --no-private --css blue --debug \
	       --navlink 'epydoc examples' doc/inh_example.py
	$(EPYDOC) -o $(HTML_EXAMPLES)/included \
	       --inheritance=included \
	       -n epydoc -u http://epydoc.sourceforge.net \
	       --no-private --css blue --debug \
	       --navlink 'epydoc examples' doc/inh_example.py
	touch .examples.up2date

# Generate the HTML version of the man page.  Note: The
# post-processing clean-up that I do is probably *not* very portable.
doc/epydoc-man.html: man/epydoc.1
	wget http://localhost/cgi-bin/man2html?epydoc -O - \
	     2>/dev/null \
	     | sed 's/<\/HEAD><BODY>/<link rel="stylesheet" href="epydoc.css" type="text\/css"\/><\/HEAD><BODY>/'\
	     | sed '/<DD>/{s/<DD>//; :loop; n; b loop;}'\
	     | sed '/<H1>/,/<HR>/{s/.*//;}'\
	     | sed 's/\(<A NAME=".*">\)&nbsp;<\/A>/\1/'\
	     | sed 's/<\/H2>/<\/H2><\/A>/'\
	     | sed 's/"\/cgi-bin\/man2html?epydocgui+1"/"epydocgui-man.html"/'\
	     | sed 's/<A HREF="\/cgi-bin\/man2html">man2html<\/A>/man2html/'\
	     > doc/epydoc-man.html

doc/epydocgui-man.html: man/epydocgui.1
	wget http://localhost/cgi-bin/man2html?epydocgui -O - \
	     2>/dev/null \
	     | sed 's/<\/HEAD><BODY>/<link rel="stylesheet" href="epydoc.css" type="text\/css"\/><\/HEAD><BODY>/'\
	     | sed '/<H1>/,/<HR>/{s/.*//;}'\
	     | sed 's/\(<A NAME=".*">\)&nbsp;<\/A>/\1/'\
	     | sed 's/<\/H2>/<\/H2><\/A>/'\
	     | sed 's/"\/cgi-bin\/man2html?epydoc+1"/"epydoc-man.html"/'\
	     | sed 's/<A HREF="\/cgi-bin\/man2html">man2html<\/A>/man2html/'\
	     > doc/epydocgui-man.html

##//////////////////////////////////////////////////////////////////////
## Standard Library docs
##//////////////////////////////////////////////////////////////////////

SLNAME = 'Python 2.2 Standard Library'
SLURL = "http://www.python.org/doc/2.2/lib/lib.html"
SLLINK = '<font size="-2">Python 2.2<br />Standard Library</font>'
SLFILES = $(shell find /usr/lib/python2.2/ -name '*.py' -o -name '*.so' \
	      |grep -v '/python2.2/config/' \
	      |grep -v '/python2.2/lib-old/' \
	      |grep -v '/python2.2/site-packages/' \
              |grep -v '/__')
export TZ='XXX00XXX;000/00,000/00' # So tzparse won't die?
stdlib-html: .stdlib-html.up2date
.stdlib-html.up2date: $(PY_SRC)
	rm -rf $(HTML_STDLIB)
	mkdir -p $(HTML_STDLIB)
	$(EPYDOC) -o $(HTML_STDLIB) -c white \
	       -n $(SLNAME) -u $(SLURL) --docformat plaintext --debug \
	       --show-imports --navlink $(SLLINK) --builtins $(SLFILES)
	touch .stdlib-html.up2date

# (this will typically cause latex to run out of resources)
stdlib-pdf: .stdlib-pdf.up2date
.stdlib-pdf.up2date: $(PY_SRC)
	rm -rf $(LATEX_STDLIB)
	mkdir -p $(LATEX_STDLIB)
	$(EPYDOC) --pdf -o $(LATEX_STDLIB) \
		--no-private -n $(SLNAME) --docformat plaintext \
		--debug --builtins $(SLFILES)
##//////////////////////////////////////////////////////////////////////
## Unit Testing
##//////////////////////////////////////////////////////////////////////

tests: test
test:
	python src/epydoc/test/__init__.py

##//////////////////////////////////////////////////////////////////////
## Other test targets
##//////////////////////////////////////////////////////////////////////

docutils: docutils-html docutils-pdf

docutils-html: .docutils-html.up2date
.docutils-html.up2date: $(PY_SRC)
	rm -rf $(HTML)/docutils
	mkdir -p $(HTML)/docutils
	$(EPYDOC) -o $(HTML)/docutils -n 'Docutils' --html \
	        --docformat plaintext --ignore-param-mismatch \
	        /usr/lib/python2.2/site-packages/docutils
	touch .docutils-html.up2date

docutils-pdf: .docutils-pdf.up2date
.docutils-pdf.up2date: $(PY_SRC)
	rm -rf $(LATEX)/docutils
	mkdir -p $(LATEX)/docutils
	$(EPYDOC) -o $(LATEX)/docutils -n 'Docutils' --pdf \
	        --docformat plaintext --ignore-param-mismatch \
	        /usr/lib/python2.2/site-packages/docutils
	touch .docutils-pdf.up2date


