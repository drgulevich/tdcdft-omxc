============
Installation
============

Download
========
Copy or clone the source from GitLab repository::

	$ git clone [the link to be provided]

Requirements
============

* C++ compiler
* `Armadillo C++ library <http://arma.sourceforge.net>`_

Installation
============

The code is build in-source with CMake. From the source folder run::

    $ cmake .

    $ make

will generate a static library ``tdcdft-omxc.a`` in the source directory.

Tests
=====
Run tests::

    $ cd tests

    $ cmake .

    $ make

    $ ./tests

Examples
========
Run examples::

    $ cd examples

    $ cmake .

    $ make

e.g. example infqw.cpp outputs dynamics of the dipole moment which can be processed as follows:

    $ ./infqw > out

    $ ./plot.py

