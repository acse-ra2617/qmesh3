
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

This directory contains unit-tests for the lib module of qmesh.

You can invoke all tests in this directory by:
```
python -m unittest discover
```

You can invoke test-cases in this directory as follows:
```
./ <filename>
```
where `<filename>` is composed as `test_<qmesh class>.py`

You can invoke an individual test routine as follows:
```
python -m unittest <test-case>.<routine>
```
where `<test-case>` is the filename without the `.py` extension.
