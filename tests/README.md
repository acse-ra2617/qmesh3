
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

This directory contains Qmesh tests.
Qmesh uses the Python `unnittest` package for testing.

There are two main layers of ontological abstraction in the Qmesh testing infrastructure:

 1. The tests, contained in a directory inside the top-most Qmesh directory. The present tests directory is a package, and thus contains an `__init__.py` file.
   This layer will be termed the *tests package* or *tests*.
 2. A classification of tests according to testing methodology and scope: unit-testing, regression-testing and so on. The categories are implemented as sub-packages.
   This layer will be termed the *categories sub-package*, or *categories*.

There are further layers and classifications within each category.
The structure of such further layers is not consistent across all categories, since the requirements, aims and therefore implementation of each category is different.
However, at the deepest levels of the overall tests structure you will always find the following three layers:

* First, a *testing module*. This is a standard python module implemented as a file that contains the next two levels: classes and functions.
* Second, a *testing class*.
* Third, the *testing functions* that are implemented as methods of the testing class.
  Note that the testing classes also define a number of members that while not directly identified as testing functions they are still useful to the testing procedure.
  Such methods/functions will set-up the environment and clean-up he workspace after testing, and therefore are run during initialisation and deletion of testing classes.

Please consult the documentation within each category directory (sub-package) for detailed information on each testing sub-package.

Contents
=========

In alphabetical order, the contents of this directory are:

 * `__init__.py`: Python initialisation file, making this directory into a package.
This file is necessary for the testing scaffolding in `setup.py` to
automatically pick-up this directory.
 * `regression`: Directory (Python package) containing Qmesh unit tests.
 * `unit`: Directory (Python package) containing Qmesh regression tests.
 * `README.md`: The present file.

Methodology
===========

You can invoke all tests in this sub-package by issuing the following command inside the sub-package directory:
```
python3 -m unittest discover
```

Individual testing modules, testing classes and testing functions can be run with commands at this directory.
Specific instructions depend on the structure of each testing sub-package.
Broadly, a file (i.e. a testing module) must be specified as an argument to the python interpreter.

Thus, when explicitly loading the `unittest` module at the command line the following template must be adopted:
```
python3 -m unittest <test package structure in dot-delimited specification>.<test module filename>
```
where `<test module filename>` is the filename without the `.py` extension.
For example, the tests in the `test_gmsh` module can be run from the present directory with:
```
python3 -m unittest unit.test_mesh.test_gmsh
```
With the dot-delimited notation you can cascade further into the testing structure and invoke an individual test class or testing function as follows:
```
python3 -m unittest <test package structure in dot-delimited specification>.<test module filename>.<test class>
python3 -m unittest <test package structure in dot-delimited specification>.<test module filename>.<test class>.<test function>
```
Where again `<test module filename>` is the filename without the `.py` extension.
 
Alternatively, tests can be run without explicitly specifying the `unittest` module at the command line, but the path to the test module must be given in standard notation.
This approach has the advantage of enabling debugging or other modules to be used during development and testing.
For example in Unix/Linux systems the command template is:
```
python3 <path to package containing the testing module>/<module file name>
```

The above example becomes:
```
python3 -m pdb unit/test_mesh/test_gmsh.py
```

Finally, note that the two approaches can be combined, with multiple modules specified at the command line:
```
python3 -m pdb -m unittest unit.test_mesh.test_gmsh
```
