from WarpVisItUtil import pError
from WarpVisItUtil import pDebug
import parallel
import multiprocessing
from visit import visit
import os

#-----------------------------------------------------------------------------
class WarpVisItCLI:
    """
    Helper class which for non-interactive runs forks a CLI
    which in turn launches a Viewer. WarpVisIt CLI manages
    the forked processes.
    """
    #-------------------------------------------------------------------------
    def __init__(self):
        """ """
        self.__CommRank = parallel.get_rank()
        self.__CLIRank = 0
        self.__CLIProc = None

    #-------------------------------------------------------------------------
    def SetSimFile(self, fileName):
        """
        Set the .sim2 file name.
        """
        self.__SimFile = fileName

    #-------------------------------------------------------------------------
    def GetSimFile(self):
        """Return the .sim2 file name"""
        return self.__SimFile

    #-------------------------------------------------------------------------
    def Initialize(self, args=['-nowin']):
        """
        Open the VisIt viewer, connect to the simulation and then execute
        the visualization defined by the visFunction. Any given arguments
        are passed on to the viewer in the form of command line arguments.
        """
        pDebug('WarpVisItCLI::Initialize')
        if self.__CommRank == self.__CLIRank:
            # launch the CLI which in turn launches the viewer.
            # it's in a separate process, this is VisIt's design.
            self.__CLIProc = multiprocessing.Process(
                  target=CLIMain,
                  args=(self.__SimFile, args))
            self.__CLIProc.daemon = True
            self.__CLIProc.start()
        return True

    def Finalize(self):
        """
        Shutdown clean up etc...
        """
        pDebug('WarpVisItCLI::Finalize')
        return

#-----------------------------------------------------------------------------
def CLIMain(simFile, args=[]):
    """
    This is the CLI main loop. The CLI resides in a separate process.
    """
    from WarpVisItUtil import pDebug
    from WarpVisItUtil import pError
    from visit import visit
    import socket

    # Were actually in the CLI, here is where we start the viewer
    # process. yes another process.
    for arg in args:
        visit.AddArgument(arg)
    visit.Launch()
    visit.OpenDatabase(simFile)

    while 1:
        pass
        # keep this process alive.

    # CLI is finished
    pDebug('CLI finished')
    exit(0)
