****************
User Interface
****************

Command Line Interface
======================

Users can run QMesh at the command line with the general command:

.. code-block:: bash

   qmesh <optional general flags> command <command-specific flags>

The optional general flags available are:

* ``-h``: Display the help message.
* ``-l``: Log info/debugging/warning/error messages to a file, specified by the user
* ``-v``: Specifies the verbosity of the logging (info, debug, warning, error, critical)

Each of the available commands and their optional flags are discussed in the following subsections

version
-------

This command displays the version number of QMesh.

git_sha_key
-----------

This command displays the revision (in the form of a SHA-1 hash) of the QMesh source code that is currently in use. This can be used, for example, when running the ``publish`` command with the ``-v`` flag.

license
-----------

This command displays the content in the ``LICENSE`` file.

generate_mesh
-------------

This command is used to create .geo and .msh files from line and polygon shape files, and raster files.

rst2fld
-------

This command is used to convert a raster file into a gmsh field.

publish
-------

This command can be used to publish the QMesh source code, as well as the input and output data files for a given project, to an online repository hosted by `Figshare <http://www.figshare.com>`_. This is useful for the purposes of recomputability and reproducibility.

The command-specific flags available are:

* ``-s``: Publish the software source code. The path to the QMesh root directory must be provided.
* ``-v``: Use this in conjunction with ``-s`` to publish a specific revision (in the form of a SHA-1 hash) of the source code. This flag is optional; without it, the revision of the source code currently in use will be published.
* ``-i``: Publish input data files (e.g. the .shp files) associated with a given project. The path to the ``.qgs`` project file must be provided.
* ``-o``: Publish output data files (e.g. the .msh files) associated with a given project. The path to the ``.qgs`` project file must be provided.
* ``-p``: Keep all publications private. This can be used in conjunction with ``-s``, ``-i`` and ``-o``. Note that the DOI associated with the publication will not be usable until it is manually made public.

For example, when publishing the QMesh software source code, use

.. code-block:: bash

   qmesh publish -s /path/to/qmesh
   
whereas for input data, use

.. code-block:: bash

   qmesh publish -i /path/to/my_project.qgs

Note that this functionality requires the `PyRDM library <https://github.com/pyrdm/pyrdm>`_ to be installed.

Upon completion of the publication process, a Digital Object Identifier (DOI) will be printed out to the screen. Such DOIs can be used to properly cite the data produced, and the precise revision of the software used to produce that data.

Graphical User Interface
========================

The graphical user interface for QMesh is in the form of a plugin for QGIS. To start the plugin, select ``QMesh`` in the ``Plugins > QMesh`` menu item. To create a Gmsh mesh file, you will need to specify the boundary vector file from the drop-down menu and the location where the mesh file should be created by clicking ``Browse files``. When ready, click the ``Generate Mesh`` button.

It is also possible to create a mesh metric raster file under ``Plugins > QMesh > Mesh metric creation``.

Python Application Programming Interface
========================================

QMesh has been split such that the main functionality can be imported into your own Python script and used programmatically. Simply import the module(s) described in the `API documentation <source/qmesh.html>`_.
