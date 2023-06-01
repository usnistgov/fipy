========
Examples
========

.. note::

   Any given module ":mod:`example.something.input`" can be found in the
   file ":file:`examples/something/input.py`".

These examples can be used in at least four ways:

- Each example can be invoked individually to demonstrate an
  application of :term:`FiPy`::

    $ python examples/something/input.py

- Each example can be invoked such that when it has finished running, you
  will be left in an interactive Python_ shell::

    $ python -i examples/something/input.py

  At this point, you can enter Python_ commands to manipulate the model or
  to make queries about the example's variable values.  For instance, the
  interactive Python sessions in the example documentation can be typed in
  directly to see that the expected results are obtained.

- Alternatively, these interactive Python_ sessions, known as doctest_
  blocks, can be invoked as automatic tests::

    $ python setup.py test --examples

  In this way, the documentation and the code are always certain to be
  consistent.

- Finally, and most importantly, the examples can be used as a templates to
  design your own problem scripts.

  .. note::

     The examples shown in this manual have been written with particular
     emphasis on serving as both documentation and as comprehensive tests
     of the FiPy framework.  As explained at the end of
     ``examples/diffusion/steadyState/mesh1D.py``, your own scripts can be
     much more succinct, if you wish, and include only the text that follows
     the "``>>>``" and "``...``" prompts shown in these examples.

     To obtain a copy of an example, containing just the script instructions,
     type::

       $ python setup.py copy_script --From x.py --To y.py


In addition to those presented in this manual, there are dozens of other
files in the :file:`examples/` directory, that demonstrate other uses of FiPy.
If these examples do not help you construct your own problem scripts,
please `contact us`_.


.. _Python:         http://www.python.org/
.. _doctest:        http://www.python.org/doc/current/lib/module-doctest.html
.. _contact us:     mailto:fipy@list.nist.gov
