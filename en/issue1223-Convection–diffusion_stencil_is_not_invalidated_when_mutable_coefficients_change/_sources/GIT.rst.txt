---------
Git usage
---------

All stages of :term:`FiPy` development are archived in a Git
repository at GitHub_. You can browse through the code at
https://github.com/usnistgov/fipy and, using a `Git client`_, you can
download various tagged revisions of :term:`FiPy` depending on your needs.

.. attention::

   Be sure to follow :ref:`INSTALLATION` to obtain all the prerequisites for
   :term:`FiPy`.

Git client
==========

A ``git`` client application is needed in order to fetch files from our
repository. This is provided on many operating systems (try executing
``which git``) but needs to be installed on many others. The sources to
build Git, as well as links to various pre-built binaries for
different platforms, can be obtained from http://git-scm.com/.

Git branches
============

In general, most users will not want to download the very latest state of
:term:`FiPy`, as these files are subject to active development and may not behave
as desired. Most users will not be interested in particular version numbers
either, but instead with the degree of code stability. Different branches are
used to indicate different stages of :term:`FiPy` development. For the
most part, we follow `a successful Git branching model`_. You will
need to decide on your own risk tolerance when deciding which stage of
development to track.

.. _cloning the repository:

A fresh copy of the :term:`FiPy` source code  can be obtained with::

   $ git clone https://github.com/usnistgov/fipy.git

An existing Git checkout of FiPy can be shifted to a different `<branch>` of
development by issuing the command::

   $ git checkout <branch>

in the base directory of the working copy. The main branches for FiPy are:

``master``
    designates the (ready to) release state of FiPy. This code is stable
    and should pass all of the tests (or should be documented that it does
    not).

Past releases of FiPy are tagged as

``x.y.z``
    Any released version of FiPy will be designated with a fixed tag: The
    current version of FiPy is |version|.  (Legacy ``version-x_y_z`` tags
    are retained for historical purposes, but won't be added to.)

Tagged releases can be found with::

   $ git tag --list

Any other branches will not generally be of interest to most users.

.. note::

   For some time now, we have done all significant development work on
   branches, only merged back to ``master`` when the tests pass
   successfully.  Although we cannot guarantee that ``master`` will never
   be broken, you can always check our :ref:`CONTINUOUSINTEGRATION` status
   to find the most recent revision that it is running acceptably.

   Historically, we merged to ``develop`` before merging to ``master``.  We
   no longer do this, although for time being, ``develop`` is kept
   synchronized with ``master``.  In a future release, we will remove the
   ``develop`` branch altogether.

For those who are interested in learning more about Git, a wide variety of
online sources are available, starting with the `official Git website`_.
The `Pro Git book`_ :cite:`ProGit` is particularly instructive.

.. _official Git website: http://git-scm.com/

.. _Pro Git book: http://git-scm.com/book

.. _GitHub: https://github.com/usnistgov/fipy

.. _a successful Git branching model: http://nvie.com/posts/a-successful-git-branching-model/
