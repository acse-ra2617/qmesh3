**********************
Setup and Installation
**********************

Dependencies
============

QMesh relies on a few other packages which must first be installed separately. These are listed in the ``README.md`` file and can be installed either using a package manager such as ``apt`` or ``pip``.

Installing from source
======================

Installation from source is currently the only way of installing qmesh, and is recommended for developers. This can be accomplished by the following steps:

* Clone the Git repository

.. code-block:: bash

   git clone https://bitbucket.org/AlexandrosAvdis/qmesh.git

* Change directory to the qmesh base directory (sometimes referred to as the qmesh `root` directory in this documentation):

.. code-block:: bash

   cd qmesh
   
* Use the Python script ``setup.py`` to install qmesh:

.. code-block:: bash

   sudo python setup.py install

Adding qmesh plugin to QGIS
-----------------------------

The first step is to deploy the QMesh plugin. Navigate to the ``qgis_plugins`` directory from the qmesh base directory using

.. code-block:: bash

   cd qgis_plugins

and type

.. code-block:: bash

   make deploy

to copy the necessary plugin files to the QGIS directory on your computer. Afterwards, load up QGIS, open the Plugins window (Plugins > Manage and Install Plugins... from the menu bar), and then select `qmesh` from the list of installed plugins.

