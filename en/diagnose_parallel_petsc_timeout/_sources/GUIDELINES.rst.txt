======================
Development Guidelines
======================

The following are the practices that have evolved over the course of
:term:`FiPy`'s development. Stray from them at your peril.

----------------------
Object Oriented Design
----------------------

Multiple inheritance is forbidden.

-----------------
Programming Style
-----------------

To the greatest extent possible, we follow :pep:`8`. We make the following
exceptions because we had already written a great deal of code before we
knew there even was a PEP 8.

 Maximum Line Length
  The preferred place to break around a binary operator is *before* the
  operator, not after it. Some examples::

      class Rectangle(Blob):

          def __init__(self, width, height,
                       color='black', emphasis=None, highlight=0):
              if (width == 0 and height == 0
                  and color == 'red' and emphasis == 'strong'
                  or highlight > 100):
                  raise ValueError("sorry, you lose")
              if width == 0 and height == 0 and (color == 'red'
                                                 or emphasis is None):
                  raise ValueError("I don't think so -- values are %s, %s"
                                   % (width, height))
              Blob.__init__(self, width, height,
                            color, emphasis, highlight)

 Package and Module Names
  Use ``mixedCase`` rather than all-``lowercase`` names.

 Function and Method Names
  Use ``mixedCase`` rather than ``lower_case_with_underscores`` names.

Any other deviations from PEP 8 should probably be corrected, but ask if
you are unsure. Even the exceptions listed above are open to discussion by
anybody willing to do the work of changing the existing code.

-------
Testing
-------

:term:`FiPy` uses :mod:`doctest` pervasively. We have not seen a case where
:mod:`unittest` would simplify things. As a heavily numerical code, care
must be taken in checking floating-point values.

Tests should not fail as a result of a missing optional install. :term:`FiPy`
defines the following :mod:`doctest` directives to specify what conditions
are needed to run the test:

  ``GMSH``
    The :program:`gmsh` executable, version 2.0 or greater, must on the
    :envvar:`PATH`.

  ``SCIPY``
    The :mod:`scipy` package must be importable.

  ``TVTK``
    The :mod:`tvtk` (or :mod:`enthought.tvtk`) package must be importable.

  ``SERIAL``
    :term:`FiPy` must be running on a single processor.

  ``PARALLEL``
    :term:`FiPy` must be running on more than one processor.

  ``PARALLEL_2``
    :term:`FiPy` must be running on two processors.

  ``PROCESSOR_0``
    Whether serial or parallel, the test will only be run on
    processor ID 0.

  ``PROCESSOR_0_OF_2``
    The test will only be run on processor ID 0 of a two processor
    parallel job.

  ``PROCESSOR_1_OF_2``
    The test will only be run on processor ID 1 of a two processor
    parallel job.

  ``PROCESSOR_0_OF_3``
    The test will only be run on processor ID 0 of a three processor
    parallel job.

  ``PROCESSOR_1_OF_3``
    The test will only be run on processor ID 1 of a three processor
    parallel job.

  ``PROCESSOR_2_OF_3``
    The test will only be run on processor ID 2 of a three processor job.

  ``LSMLIB``
    The :mod:`lsmlib` package must be importable.


Further directives can be defined using
:func:`~fipy.tests.doctestPlus.register_skipper`. These definitions should
be placed in the file to which they pertain, e.g., ``GMSH`` is defined in
:file:`fipy/meshes/gmshImport.py`.
