##########################################
Testing standard, biblist and missing refs
##########################################

First [Brett2002]_.

Looking for [Aston2006]_ and [Brett2001]_.  See also :ref:`test2`.

.. biblisted:: some.bib
    :sort:
    :encoding: utf-8

    Brett2001
    Bishop2004

Lorem ipsum [Ref]_ dolor sit amet.

Some text

.. [Ref] Book or article reference, URL or whatever.

More text

.. [UnusedRef] A reference to which there is no reference

More text again

.. .. [ref_elsewhere] A reference from another document

A

lot

of

text

lines

to

show

that

the

reference

links

are

working

.. bibmissing:: other.bib, some.bib
    :style: test-style
    :sort:
