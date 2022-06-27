#!/bin/sh

#    Copyright (C) 2013 Alexandros Avdis and others. See the AUTHORS file for a full list of copyright holders.
#
#    This file is part of QMesh.
#
#    QMesh is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    QMesh is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with QMesh.  If not, see <http://www.gnu.org/licenses/>.

.PHONY: qmesh doc test uninstall

input: 	clean doc

install: qmesh

qmesh:
	@echo ">>> Installing qmesh..."
	pip install .

doc:
	@echo ">>> Building the documentation..."
	cd doc && make html && cd ..
	@echo "Done!"

test:
	@echo ">>> Running tests..."
	cd tests && make && cd ..
	@echo "Done!"

clean:
	@echo ">>> Cleaning up..."
	cd doc && make clean && cd ..
	cd tests && make clean && cd ..
	python setup.py clean
	@echo "Done!"

devinstall:
	@echo ">>> Installing qmesh in 'developer mode'..."
	pip install -e .

uninstall:
	@echo ">>> Un-installing qmesh..."
	pip uninstall qmesh
