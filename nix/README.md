# FiPy Nix Installation

This directory contains a Nix expression for installing FiPy using
Nix. [Nix is a powerful package manager for Linux and other Unix
systems that makes package management reliable and
reproducible](https://nixos.org/nix/).

## Getting Started with Nix

There are a number of tutorials on getting started Nix. The page that
I used when getting started is on the Blog of the HPC team of GRICAD,

 - https://gricad.github.io/calcul/nix/tuto/2017/07/04/nix-tutorial.html

I also made my own notes,

 - https://github.com/wd15/nixes/blob/master/NIX-NOTES.md

which are a succinct steps that I use when setting up a new system with
Nix.

## Installing FiPy

Once you have a working Nix installation use

    $ nix-shell --pure

in the base FiPy directory to install FiPy with Python 2. To install with
Python 3, use

    $ nix-shell --pure nix/py3.nix

## Additional Packages

To install additional packages available from Nixpkgs include them in
the `buildInputs` list in (`nix/build.nix`)[nix/build.nix].

## Using Pip

Packages unavailable from Nix can be installed using Pip. Uncomment
the `postShellHook` section in `nix/build.nix` to enable this. In this
case, the installation has been set up so that the shell knows about a
`.local` directory in the FiPy base directory used by `pip install
--user package` for installation.  So, for example, to install the
toolz package from within the Nix shell use

    $ pip install --user toolz

This is an impure action and will persist after the shell has been
closed so the `.local` directory should be deleted after use.