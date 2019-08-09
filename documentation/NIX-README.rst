Nix Installation
================

:term:`FiPy` now has a `Nix`_ expression for installing :term:`FiPy`
using `Nix`_. `Nix`_ is a powerful package manager for Linux and other
Unix systems that makes package management reliable and
reproducible. The recipe works on both Linux and Mac OS X.

Getting Started with Nix
------------------------

There are a number of tutorials on getting started with `Nix`_. The
page that I used when getting started is on the Blog of the HPC team
of GRICAD,

https://gricad.github.io/calcul/nix/tuto/2017/07/04/nix-tutorial.html

I also made my own notes,

https://github.com/wd15/nixes/blob/master/NIX-NOTES.md

which are a succinct steps that I use when setting up a new system with
Nix.

Installing
----------

Once you have a working Nix installation use::

    $ nix-shell --pure

in the base :term:`FiPy` directory to install :term:`FiPy` with Python
3 by default. Modify the `shell.nix` file to use another version of
Python. ``nix-shell`` drops the user into a shell with a working
version of :term:`FiPy`. To test your installation use::

    $ nix-shell --pure --command "python setup.py test"

.. note::

   :term:`Trilinos` is currently not available as part of the Nix
   :term:`FiPy` installation.


Additional Packages
-------------------

To install additional packages available from Nixpkgs_ include them in
the `nativeBuildInputs` list in `shell.nix`.


Using Pip
---------

Packages unavailable from Nix can be installed using :term:`Pip`. In
this case, the installation has been set up so that the Nix shell
knows about a ``.local`` directory in the base :term:`FiPy` directory
used by :term:`Pip` for installation.  So, for example, to install the
``toolz`` package from within the Nix shell use::

    $ pip install --user toolz

The ``.local`` directory will persist after the Nix shell has been
closed.

.. _Nix: https://nixos.org/nix/
.. _Nixpkgs:  https://nixos.org/nixpkgs/
