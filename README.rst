qmesh3 
=======
Qmesh3 has been created as a temporary solution to allow qmesh to function with python3 and QGIS3 within a functional pip install to facilitate users access until the original qmesh updates are released. This repository has been created by `Raul Adriaensen <https://www.linkedin.com/in/rauladriaensen/>`_ , contact email can be found `here <https://www.imperial.ac.uk/people/raul.adriaensen17>`_

qmesh is a software package for creating high-quality meshes using `QGIS <https://www.qgis.org>`_ and `Gmsh <https://geuz.org/gmsh>`_.
The meshes can be used in finite element numerical models such as `TELEMAC <http://www.opentelemac.org>`_, `Fluidity <https://www.fluidity-project.org>`_ and `Thetis <https://thetisproject.org/>`_.
For more information please visit the project `website <https://www.qmesh.org>`_.



Development, Maintenance and Licence
------------------------------------

qmesh was developed and is maintained by `Alexandros Avdis <https://orcid.org/0000-0002-2695-3358>`_ and `Jon Hill  <https://orcid.org/0000-0003-1340-4373>`_.
Please see file `AUTHORS.md <https://bitbucket.org/qmesh-developers/qmesh-containers/raw/HEAD/AUTHORS.md>`_ for more information.

qmesh is an open-source project distributed under the GNU General Public Licence v3 (`GPL v3 <https://www.gnu.org/licenses/gpl-3.0.en.html>`_).
The file `LICENCE <https://bitbucket.org/qmesh-developers/qmesh-containers/raw/HEAD/LICENSE>`_ contains a complete copy of the licence.

The source code is freely available to download and adapt, under the conditions stated in the `LICENCE <https://bitbucket.org/qmesh-developers/qmesh-containers/raw/HEAD/LICENSE>`_ file.
However, we follow a closed-development approach, where maintenance and development is carried out by the package originators.
We have opted for this approach to limit the resources needed for development: A larger development team necessitates larger management structures and significant CI/CD expenditure.
Larger management structures will absorb time intended for other exciting and useful research.
CI/CD expenditure could threaten our aim to offer a focused package, at no monetary costs to the users.



Documentation 
---------------

You can access relevant documentation through the following avenues:

* The `qmesh website <https://www.qmesh.org>`_.
* The `qmesh synoptic manual <https://qmesh-synoptic-manual.readthedocs.io/en/latest>`_.

Walk Through
---------------

*Installation*:

1. install `GD-basisChangeTools/ <https://pypi.org/project/GFD-basisChangeTools/>`_
2. install `gmsh <https://installati.one/ubuntu/20.04/gmsh/>`_
3. install `qgis <https://qgis.org/en/site/forusers/alldownloads.html>`_, however some people seem to have more success first doing "sudo apt update"\" and then "Sudo apt install qgis"
4. install `qmesh3 <https://pypi.org/project/qmesh3/>`_
5. testing installation (optional) > "python -m unittest discover <test_directory>"

*tutorials*:

All key data files and videos have been placed together under this `Onedrive <https://1drv.ms/u/s!AglgFElvf_OWl8gIx0FxAIcdOhUv8g?e=VrIak0>`_. . A brief description of the files:

-   `Semi automated Mesh Generation of the Coastal Oceans with Engineering Structures Using Satellite Data <https://www.dropbox.com/s/1bwrwjl51cnhhju/Semi-automated%20Mesh%20Generation%20of%20the%20Coastal%20Oceans%20with%20Engineering%20Structures%20Using%20Satellite%20Data.pdf?dl=0>`_
   - great tutorial to get started, this was originally run on docker but can now be done without if the above installation is followed
- additioanl_tutorials_and_docker_installation
   - this material helped me to get going and explore all the options from the docker, qmesh-cli and learn some qgis tricks
- installating_qmesh_manually
   - in case the pip install breaks, instructions on how to setup qmesh manually can be found here
- tutorials_Angeloudis_Athnasios
   - an additional set of tutorials shared 
  
Please read the relevant readme files in the OneDrive and the email chains to find the full information and names of those who kindly put all this content together

