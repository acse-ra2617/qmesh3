
[comment]: # (Copyright C 2013 Alexandros Avdis and others.)
[comment]: # (See the AUTHORS file for a full list of copyright holders.)
[comment]: # ( )
[comment]: # (This file is part of Qmesh.)
[comment]: # ( )
[comment]: # (Qmesh is free software: you can redistribute it and/or modify)
[comment]: # (it under the terms of the GNU General Public License as published by)
[comment]: # (the Free Software Foundation, either version 3 of the License, or)
[comment]: # ([at your option] any later version.)
[comment]: # ( )
[comment]: # (Qmesh is distributed in the hope that it will be useful,)
[comment]: # (but WITHOUT ANY WARRANTY; without even the implied warranty of)
[comment]: # (MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the)
[comment]: # (GNU General Public License for more details.)
[comment]: # ( )
[comment]: # (You should have received a copy of the GNU General Public License)
[comment]: # (along with Qmesh.  If not, see <http://www.gnu.org/licenses/>.)

[comment]: # (The rendering of parts of this file can be tested with https://johnmacfarlane.net/babelmark2/)
[comment]: # (Please place each sentence on a separate line, the renderer will construct the paragraphs.)

Overview
==========

This directory contains unit-tests for Qmesh.
Following the rationale presented in the documentation of the tests, the unit-tests are structured as six-layers:

 1. The *tests*, contained in a directory inside the top-most qmesh directory.
    The tests directory is a package, and thus contains an `__init__.py` file.
 2. A classification of tests according to testing methodology and scope: unit-testing, regression-testing and so on.
    The *categories* are implemented as sub-packages.
    The present category is the *unit-testing sub-package*.
 3. A further categorisation is made within each unit-testing sub-package, where each category is aimed at testing one of the main qmesh packages: vector, raster, mesh and publish.
    Each category is implemented as a separate sub-sub-package, in a separate directory.
    The directories follow the naming convention `test_<Qmesh package>`.
    The present sub-sub-package is aimed to test the `raster` qmesh package.
 4. *Testing modules*: Inside each sub-sub-package a python file includes the next three levels of the Qmesh unit-tests organisation, starting with testing modules.
    Each testing module typically tests the corresponding Qmesh module.
    The file names follow the convention `test_<Qmesh module>.py`
 5. *Testing classes*: Similar to levels above, testing classes typically test the corresponding qmesh class.
 6. *Testing functions*: Finally, testing functions typically test the corresponding Qmesh function/method.

See the `README.md` file in the `tests` directory for more information.

Contents
=========

In alphabetical order, the contents of this directory are:

 * `__init__.py`: Python initialisation file, making this directory into a package.
This file is necessary for the testing scaffolding in `setup.py` to
automatically pick-up this directory.
 * `test_raster.py`: Python module containing unit tests of the Qmesh `raster` module.
 * `README.md`: The present file.

Methodology
===========

You can invoke all tests in this sub-package by issuing the following command inside the sub-package directory:
```
python3 -m unittest discover
```

The python debugger can be used to run tests:
```
python3 -m gdb <test module filename>
```

See the `README.md` file in the `tests` directory for more information, including how to run individual testing-classes and testing-functions. 
However, keep in mind that the paths and module-notation must be relative to the directory where the command is issued.
