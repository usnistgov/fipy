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

To install additional packages with Nix add them to the `buildInputs`
list in (`nix/build.nix`)[nix/build.nix].

## Using Pip

The install has been set up to work with pip. The installation has
been set up with a `.local` directory in the base FiPy directory. So,
when using pip, run it from the base fipy directory or specify the
full path name. For example,

    $ echo $PWD
    /base/fipy/directory
    $ pip install --install-option="--prefix=$PWD/.local" my_package

will install `my_package` in `.local`. There are alternate way to set
up Pip with Nix including use virtualenv or the `--user` option. However
if you use the `--user` option or virtualenv outside of