Nix Installation
================

:term:`FiPy` now has a `Nix`_ expression for installing :term:`FiPy`
using `Nix`_. `Nix`_ is a powerful package manager for Linux and other
Unix systems that makes package management reliable and
reproducible. The recipe works on both Linux and Mac OS X. Go to
`nix.dev`_ to get started with Nix.

Installing
----------

Once you have a working Nix installation use::

    $ nix develop

in the base :term:`FiPy` directory to install :term:`FiPy` with Python
3 by default. ``nix develop`` drops the user into a shell with a working
version of :term:`FiPy`. To test your installation use::

    $ nix develop --command bash -c "python setup.py test"

.. note::

   The SciPy solvers are the only available solvers currently.


.. _Nix: https://nixos.org/nix/
.. _Nixpkgs:  https://nixos.org/nixpkgs/
.. _nix.dev: https://nix.dev
