#!/usr/bin/env python

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


import logging
LOG = logging.getLogger(__package__)


class Publisher:

   def __init__(self, path, service="figshare"):
      """ Initialise the qmesh publishing tool.
      
      :param str path: If the software is being published, this is the path to the directory
                       containing the source code. If a dataset is being published, this is
                       the path to the .qgs project file.
      :returns: None
      """
      LOG.debug("Initialising qmesh publishing tool...")
      #RDM is currently an optional. try to import PyRDM and raise exception if not found
      try:
         from pyrdm.publisher import Publisher as pyrdmPublisher
      except ImportError:
         msg = "The PyRDM Publisher class could not be imported. Research Data Management" + \
               " is an optional qmesh feature, requiring PyRDM. Please install PyRDM separately."
         LOG.error(msg, exc_info=True)
      raise ImportError(msg)
      self.publisher = pyrdmPublisher(service=service)
      #Set the path.
      self.path = path
      LOG.debug("qmesh publishing tool initialised.")
      return

   def publish(self, publication_type, version=None, private=True):
      """ Publish the qmesh source code or any project datasets.
      
      :param str publication_type: A single character representing the type of research output to be published (software, input data, or output data)
      :param str service: The online repository service that should be used (e.g. Figshare or Zenodo)
      :param str version: If it is the software that is being considered, then this is the version of the source code that should be published, in the form of a SHA-1 hash from the Git version control system. By default, this is set to None, such that the HEAD revision will be used.
      :param bool private: If True, the publication will remain private. By default, this is set to False such that the publication will be made public.  
      
      :rtype: tuple
      :returns: A tuple containing the publication ID and the corresponding DOI.
      """
      import os
      import sys
      from xml.dom import minidom
      # Publish the software
      if(publication_type == "s"):
         LOG.info("Publishing the qmesh source code...")
       
         # Get the SHA-1 hash of the software version.
         if(version is not None):
            LOG.info("Using the software version provided: %s" % version)

         # The 'pid' is the publication ID (e.g. on Figshare) and 'doi' is the Digital Object Identifier.
         pid, doi = self.publisher.publish_software(name="qmesh", local_repo_location=self.path, version=version, private=private)

      else:
         # Publish the input/output data
         LOG.info("Publishing data...")
      
         # Parse the XML
         xml = minidom.parse(self.path)
         
         # Get the root element which contains the QGIS project name.
         root = xml.getElementsByTagName('qgis')[0]
       
         project = root.attributes["projectname"].value
         if project == "": #Qgis can give an empty project-name, but figshare does not like it when used as a tag
            project = self.path
         LOG.info("Project name: %s" % project)
       
         tags = ["qmesh", project]
         if(publication_type == "i"):
            temp = "input"
            tags.append("input data")
         elif(publication_type == "o"):
            temp = "output"
            tags.append("output data")
         title = "qmesh %s data for QGIS project %s" % (temp, project)
      
         files = []
         
         if(publication_type == "i"):
            # Get all the (input) data sources from the QGIS project file.
            for element in xml.getElementsByTagName('datasource'):
               value = str(element.childNodes[0].nodeValue)
               if(os.path.dirname(self.path) != ""):
                  value = os.path.dirname(self.path) + "/" + value
               files.append(value)
         
            LOG.info("Current list of input files to publish: %s" % files)

            # Add in any additional files that the user wants to publish
            extras = raw_input("""
Please enter any additional (relative) paths to input files that you want to publish.
The paths should be provided in a Python list of strings (e.g. [\"file1.txt\", \"file2.txt\"]).
If there are no additional files, just press Enter to continue.
""")
            if(extras != ""):
               files = files + eval(extras)
               
         elif(publication_type == "o"):
            # Get all output data file paths directly from the user.
            files = raw_input("""
Please enter all (relative) paths to the output files that you want to publish.
The paths should be provided in a Python list of strings (e.g. [\"file1.txt\", \"file2.txt\"]).
If there are no additional files, just press Enter to continue.
""")
            if(files != ""):
               files = eval(files)
            else:
               LOG.error("No output files specified.")
               sys.exit(1)

         LOG.info("Files to publish: %s" % files)
         
         # Publish the datas.
         parameters = {"title":title, "description":title, "files":files, "category":"Computational Physics", "tag_name":tags}
         pid, doi = self.publisher.publish_data(parameters, pid=None, private=private)
         
      LOG.info("Finished publishing.")
      return pid, doi

