--------
Examples
--------

These examples can be used in at least four ways:

- Each example can be invoked individually to demonstrate an 
  application of |FiPy|::

    $ examples/something/input.py

- Each example can be invoked such that when it has finished running, you
  will be left in an interactive Python_ shell::

    $ python -i examples/something/input.py

  At this point, you can make Python_ queries about the example's variable
  values.  For instance, the interactive Python sessions in the example
  documentation can be typed in directly to see that the expected results
  are obtained.

- Alternatively, these interactive Python_ sessions, known as doctest
  blocks, can be invoked as automatic tests::

    $ python setup.py test --examples

  In this way, the documentation and the code are always certain to be
  consistent.

- Finally, the examples can be used as a templates to design your own problem
  scripts.


.. _Python:               http://www.python.org/
