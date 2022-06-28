#    Copyright (C) 2017 Alexandros Avdis and others. See the AUTHORS file for a full list of copyright holders.
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


# Logging-related functions to set up a logger which outputs coloured text. This is based on the code by 'airmind' on Stack Overflow:
# http://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output

import logging
import qgis.core

_BLACK, _RED, _GREEN, _YELLOW, _BLUE, _MAGENTA, _CYAN, _WHITE = range(8)

#These are the sequences need to get colored ouput
_RESET_SEQ = "\033[0m"
_COLOR_SEQ = "\033[1;%dm"
_BOLD_SEQ = "\033[1m"

_COLORS = {
    'WARNING': _YELLOW,
    'INFO': _WHITE,
    'DEBUG': _BLUE,
    'CRITICAL': _MAGENTA,
    'ERROR': _RED
}

def formatter_message(message, use_color = True):
    r""" Return a formatted version of the message.
    
    :param str message: The non-formatted message.
    :param bool use_color: Use a colour format for the message.
    :returns: The formatted version of the message.
    :rtype: str
    """
    if use_color:
        message = message.replace("$RESET", _RESET_SEQ).replace("$BOLD", _BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message

class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, use_color = True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in _COLORS:
            levelname_color = _COLOR_SEQ % (30 + _COLORS[levelname]) + levelname + _RESET_SEQ
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)
        
# Custom logger with multiple destinations
def ColoredLogger(name):
    FORMAT = "[$BOLD%(name)s$RESET] %(levelname)s: %(message)s"
    COLOR_FORMAT = formatter_message(FORMAT, True)
    logger = logging.getLogger(name)                
    logger.setLevel(logging.DEBUG)
    color_formatter = ColoredFormatter(COLOR_FORMAT)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(color_formatter)
    logger.addHandler(consoleHandler)
    return logger

LOG = ColoredLogger("qmesh")
      
def setLogOutputFile(path):
    FORMAT = "[$BOLD%(name)s$RESET] %(levelname)s: %(message)s"
    COLOR_FORMAT = formatter_message(FORMAT, True)
    color_formatter = ColoredFormatter(COLOR_FORMAT)
    file_handler = logging.FileHandler(path)
    file_handler.setFormatter(color_formatter)
    LOG.addHandler(file_handler)
    LOG.info('qmesh log output also copied to file '+path)
      
class _subprocess_log_queue():
  def __init__(self, subproc, subproc_name, stdoutFileName, stderrFileName):
    import queue
    import threading
    self.stdOutErrQueue = queue.Queue()
    self.subproc = subproc
    self.stdoutFileName = stdoutFileName
    self.stderrFileName = stderrFileName
    stdout_parser_thread = threading.Thread(target=self.stdout_parser,
                                            name=subproc_name+'_stdout_parser')
    stderr_parser_thread = threading.Thread(target=self.stderr_parser,
                                            name=subproc_name+'_stderr_parser')
    stdout_parser_thread.start()
    stderr_parser_thread.start()
  def stdout_parser(self):
    '''Capture subprocess stdout and store in queue object'''
    stdoutFile = open(self.stdoutFileName,'r')
    while True:
      line = stdoutFile.readline()
      if len(line.split()) != 0:
        self.stdOutErrQueue.put(line)
      elif len(line.split()) == 0 and self.subproc.poll() is not None:
        stdoutFile.close()
        break
  def stderr_parser(self):
    '''Capture subprocess stderr and store in queue object'''
    stderrFile = open(self.stderrFileName, 'r')
    while True:
      line = stderrFile.readline()
      if len(line.split()) != 0:
        self.stdOutErrQueue.put(line)
      elif len(line.split()) == 0 and self.subproc.poll() is not None:
        stderrFile.close()
        break
  def emptyQueue(self):
     '''Check the queue storing stdout and stderr messages, and return True if it is empty. Return False if not empty'''
     return self.stdOutErrQueue.empty()
  def getLine(self):
     line = self.stdOutErrQueue.get(True,0.5)
     return line

def initialise():
    LOG.warning("The initialise method is deprecated and will be removed in a future version.")
    LOG.warning("Initialisation is now carried out when qmesh is imported.")
    LOG.warning("Explicit calls of the initialise method by the users are no longer required.")
    LOG.warning("You can delete calls of the initialise method in user scripts.")
    LOG.warning("For provenance information use qmesh.provenance()")

def provenance():
    """ List provenance information.

    List provenance information (version, and git-SHA-key) for qmesh, as well
    as for the main dependencies: qgis (version and initialisation files location)
    and gmsh (version and location of CLI binary).
    """
    from .__init__ import __version__
    from .__init__ import __git_sha_key__
    from .__init__ import __qgis_path__
    from .__init__ import __gmsh_bin_path__
    #
    LOG.info('Qmesh version ' + __version__ + ' with git sha key ' + __git_sha_key__)
    if __version__ is None or __git_sha_key__ is None:
        LOG.info('qmesh probably imported from non-packaged distribution.')
    # Look for qgis and report version and location
    # From qgis2 to qgis3 the location of QGIS_VERSION changed
    try:
        qgis_version = qgis.core.QGis.QGIS_VERSION
    except AttributeError:
        qgis_version = qgis.core.Qgis.QGIS_VERSION
    LOG.info('QGIS version ' + qgis_version + ' installed at ' + __qgis_path__)
    # Look for gmsh and report version and location
    LOG.info('gmsh installed at ' + __gmsh_bin_path__)
