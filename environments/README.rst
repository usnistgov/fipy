.. _CREATE_CONDA_ENVIRONMENT:

Create a conda_ environment
===========================

Use one of the following methods to create a self-contained conda_
environment and then download and populate the environment with the
prerequisites for :term:`FiPy` from the conda-forge_ channel at
https://anaconda.org.  See `this discussion
<https://pythonspeed.com/articles/conda-dependency-management/>`_
of the merits of and relationship between the different methods.

* Conda_ environment files

  This option is the most upgradable in the future and probably the best
  for development.

  ::

    $ conda env create --name <MYFIPYENV> \
        --file environments/<SOLVER>-environment.yml

  .. note::

     You can try to include multiple solver suites using ``conda env
     update``, but be aware that different suites may have incompatible
     requirements, or may restrict installation to obsolete versions of
     Python.  Given that :term:`FiPy` can only use one solver suite during
     a run, installing more than one solver in an environment isn't
     necessary.

     .. attention::

        Successively updating an environment can be unpredictable, as later
        packages may conflict with earlier ones.  Unfortunately, ``conda
        env create`` `does not support multiple environment files
        <https://github.com/conda/conda/issues/9294>`_.

        Alternatively, combine the different
        :file:`environments/<SOLVER>-environment.yml` files you wish to
        use, along with `environment.yml` files for any other packages you
        are interested in (`conda-merge
        <https://github.com/amitbeka/conda-merge>`_ may prove useful).
        Then execute::

          $ conda env create --name <MYFIPYENV> --file <MYMERGEDENVIRONMENT>.yml

* conda-lock_ lockfiles

  This option will pin all the packages, so is the most reproducible, but
  not particularly upgradable.  For most, this is the safest way to
  generate a FiPy environment that consistently works.

  ::

    $ conda-lock install --name <MYFIPYENV> \
        environments/locks/conda-<SOLVER>-lock.yml

  or, to be really explicit (and obviating the need for conda-lock_)::

    $ conda create --name <MYFIPYENV> \
        --file environments/locks/conda-<SOLVER>-<PLATFORM>.lock

* Directly from conda-forge_, picking and choosing desired packages

  This option is the most flexible, but has the highest risk of missing or
  incompatible packages.

  e.g.::

    $ conda create --name <MYFIPYENV> --channel conda-forge \
        python=3 numpy scipy matplotlib-base future packaging mpich \
        mpi4py petsc4py mayavi "gmsh <4.0|>=4.5.2"

  or::

    $ conda create --name <MYFIPYENV> --channel conda-forge \
        python=2.7 numpy scipy matplotlib-base future packaging \
        pysparse mayavi "traitsui<7.0.0" "gmsh<4.0"

  .. attention::

     Bit rot has started to set in for Python 2.7.  One consequence is that
     :class:`~fipy.viewers.vtkViewer.VTKViewer`\s can raise errors
     (probably other uses of :term:`Mayavi`, too). Hence, the constraint
     of `"traitsui<7.0.0"`.

.. _conda-lock: https://github.com/conda/conda-lock

