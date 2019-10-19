# tdcdft-omxc

See detailed documentation at GitLab pages: https://drgulevich.gitlab.io/tdcdft-omxc


**Requirements**

    C++ compiler

    Armadillo C++ library

**Installation**

The code is build in-source with CMake. From the source folder run:

    $ cmake .
    $ make

will generate a static library tdcdft-omxc.a in the source directory.

**Tests**

Run tests:

    $ cd tests
    $ cmake .
    $ make
    $ ./tests

**Examples**

Run examples:

    $ cd examples
    $ cmake .
    $ make

e.g. example qw.cpp outputs dynamics of the dipole moment which can be processed as follows:

    $ ./qw > out
    $ ./plot.py

**Publications**

1. D. R. Gulevich, Ya. V. Zhumagulov, A. V. Vagov, and V. Perebeinos, Nonadiabatic time-dependent density-functional theory at the cost of adiabatic local density approximation, https://arxiv.org/abs/1908.10941.

