import os
import sys
import time
import argparse
from WarpVisItUtil import pError,pDebug,pStatus,getEnvVar

#-----------------------------------------------------------------------------
class WarpVisItCLI:
    """
    Helper class which for non-interactive runs forks a CLI
    which in turn launches a Viewer. WarpVisIt CLI manages
    the forked processes.
    """
    #-------------------------------------------------------------------------
    def __init__(self, args=[]):
        """ """
        self.__Timeout = 1000
        self.__Interval = 10
        self.__SimFile =  getEnvVar('WARPVISIT_SIM2_FILE',str,'WarpVisIt.sim2')
        self.__Log = getEnvVar('WARPVISIT_CLI_LOG',str,'')
        self.__ViewerOpts = None

        # parse command line args
        ap = argparse.ArgumentParser(usage=argparse.SUPPRESS,prog='WarpVisItCLI',add_help=False)
        ap.add_argument('--sim-file',type=str,default=self.__SimFile)
        ap.add_argument('--timeout',type=int,default=self.__Timeout)
        ap.add_argument('--cli-log',type=str,default=self.__Log)
        ap.add_argument('--viewer-opts',type=str,default=None)
        opts = vars(ap.parse_known_args(args)[0])
        self.__SimFile = os.path.abspath(opts['sim_file'])
        self.__Timeout = opts['timeout']
        self.__Log = opts['cli_log']
        if self.__Log:
            f = open(self.__Log,'w')
            sys.stderr = f
            sys.stdout = f
        self.__ViewerOpts = opts['viewer_opts']
        return

    #-------------------------------------------------------------------------
    @staticmethod
    def GetCommandLineHelp():
        """
        Return list of command line options and environment vars
        """
        return ("--sim-file : WARPVISIT_SIM2_FILE : Path to read sim2 file\n"
          "--timeout : : Number of seconds to wait for the engine\n"
          "--cli-log : : Redirect stderr and stdout to this file\n"
          "--viewer-opts : : Command linie options to pass to the VisIt viewer\n")


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
    def Initialize(self):
        """
        Configure the object. If the sim2 file is
        located before timeout occurs True is returned.
        """
        if __debug__: pDebug('WarpVisItCLI::Initialize')
        # we may need to wait while the engine launches
        # in a separate process
        simFileFound=False
        n = self.__Timeout/10
        i = 0
        while i<n:
            if os.path.isfile(self.__SimFile):
                simFileFound=True
                sys.stderr.write('CLI found sim file (%s)'%(self.__SimFile))
                break
            else:
                sys.stderr.write('.')
                time.sleep(10)
            i += 1
        sys.stderr.write('\n')

        if not simFileFound:
            pError('CLI failed to find sim file (%s)'%(self.__SimFile))
            return False

        return True

    #-------------------------------------------------------------------------
    def EventLoop(self, args=['-nowin']):
        """
        Open the VisIt viewer, connect to the simulation and then execute
        the visualization defined by the visFunction. Any given arguments
        are passed on to the viewer in the form of command line arguments.
        """
        from WarpVisItUtil import pDebug,pError
        from visit import visit
        import sys
        import os
        if __debug__: pDebug('WarpVisItCLI::EventLoop')

        # this process becomes the CLI. here is where we start the viewer
        # process.
        if self.__ViewerOpts:
            args += self.__ViewerOpts.split()

        for arg in args:
            visit.AddArgument(arg)

        visit.Launch()
        ok = visit.OpenDatabase(self.__SimFile)

        # rm the sim file. if left hanging around
        # it could cause the next run to fail and
        # we don't need it after we've opened it
        if ok:
            os.unlink(self.__SimFile)

        while ok:
            time.sleep(500)
            pass
            # keep this process alive.

        # CLI is finished
        pError('CLI failed to launch')
        return False

    #-------------------------------------------------------------------------
    def Finalize(self):
        """
        Shutdown clean up etc...
        """
        sys.stderr.write('WarpVisItCLI normal termination\n')
        return
