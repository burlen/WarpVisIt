import os
import sys
import argparse
import warp
from WarpVisItUtil import pError
from WarpVisItUtil import pDebug

#############################################################################
def NewWarpVisItSimulation(args=[]):
    """
    example of the user supplied factory method that WarpVisIt uses
    to create an object instance that implements WarVisItSimulation
    interface. args are command line arguments perhaps as obtained
    from sys.argv
    """
    return AWarpVisItSimulationImplementaion(args)

#############################################################################
class WarpVisItSimulation(object):
    """
    WarpVisItSimulation

    Class describing interface used by WarpVisIt to control a running
    simulation. This class provides a default implementation for most
    of the methods. At a minimum a user will need to override Initialize
    to handle run specific warp iniziation. Users expose their implementation
    to WarpVisIt through a factory function which should defined at
    file scope in the file containing their WarpVisItSimulation implementation
    and be named NewWarpVisItSimulation(args), where args is a list of command
    line arguments (could be sys.argv). The run time is made aware of these
    through the --warp-script command line option, which should contain
    the path to the file. For example a user want to run a simulation called
    MySim, then they must at a minimum provide a file with the following:

    ##########################
    ##### begin MySim.py #####

    # 1: define the factory
    def NewWarpVisItSimulation(args=[])
        return MySim(args)

    # 2: implement the simulation class interface
    class MySim(WarpVisItSimulation):
        def __init__(self, args=[]):
            WarpVisItSimulation.__init__(self,args)
            self.AddRenderScript('plotV','plotV.py')
            # etc ...

        def Initialize(self):
            # setup and start warp here
            WarpVisItSimulation.Initialize(self)

        # override other methods if needed.

    ##### end MySim.py #####
    ########################

    and start WarpVisIt with rthe command line options:
    --warp-script=/some/path/MySim.py --script-dir=/path/to/dir/of/poltV.py
    """

    #-------------------------------------------------------------------------
    def __init__(self, args=[]):
        """
        Initialize the object and optionaly process common command
        line arguments.
        """
        self._ScriptDir = os.getenv('WARPVISIT_SCRIPT_DIR')
        if not self._ScriptDir:
            self._ScriptDir = ''
        self._Start = 0
        self._Stop = 500
        self._Step = 100
        self._Plot = 100
        self._RenderScripts = None
        self._HaveRenderScripts = False

        # parse command line args
        ap = argparse.ArgumentParser(usage=argparse.SUPPRESS,prog='WarpVisItSimulation',add_help=False)
        ap.add_argument('--script-dir',type=str,default=self._ScriptDir)
        ap.add_argument('--start',type=int,default=self._Start)
        ap.add_argument('--stop',type=int,default=self._Stop)
        ap.add_argument('--step',type=int,default=self._Step)
        ap.add_argument('--plot',type=int,default=self._Plot)
        opts = vars(ap.parse_known_args(args)[0])
        self._ScriptDir = os.path.abspath(opts['script_dir'])
        self._Start = opts['start']
        self._Stop = opts['stop']
        self._Step = opts['step']
        self._Plot = opts['plot']

        # handle special case for no plots.
        if self._Plot == 0:
            self._RenderScripts = {}
            self._HaveRenderScripts = True

        return

    #-------------------------------------------------------------------------
    @staticmethod
    def GetCommandLineHelp():
        """
        Return string listing command line options and environment
        vars(if any) one command per line. This is used in repsonse
        to --help
        """
        return ("--script-dir : : Path to directory containing vis scripts\n"
            "--start : : Number of steps to take before visualization\n"
            "--stop  : : Step number to end the run at\n"
            "--step  : : Number of steps to advance with each update\n"
            "--plot  : : Number of steps between visualizations\n")

    #-------------------------------------------------------------------------
    def GetRenderScripts(self):
        """
        This function returns a dictionary of rendering scripts
        whose key is a descriptive string. This dictionary will
        be used when our Render function returns a list of scripts
        to run to access the script source code. The source is
        given to VisIt CLI for execution.

        This implementation expects that the user will have called
        AddRenderScripts for each script, initializing the key and
        script scource file name. The first time we are called we
        load the scripts replacing source file with the source
        itself.
        """
        # return the disctionary
        if self._HaveRenderScripts:
            return self._RenderScripts

        # transform the dictionary values from source file names
        # into source code.
        if self._RenderScripts is None:
            pError('RenderScripts uninitialized')
            return {}

        for key,fileName in self._RenderScripts.iteritems():
            if not self._ScriptDir is None:
                fileName = os.path.join(self._ScriptDir,fileName)
            f = open(fileName)
            code = f.read()
            f.close()
            self._RenderScripts[key] = code

        self._HaveRenderScripts = True

        return self._RenderScripts

    #------------------------------------------------------------------------
    def GetActiveRenderScripts(self):
        """
        If the current time step should be visualized return a list
        of keys naming which rendering scripts should run. If the
        list is empty then nothing will be rendered this time step.

        This implementation will return all keys in the dictionary
        produced by GetRenderScripts every n steps, where n can be
        passed on the command line via --plot option or by calling
        SetPlotInterval.
        """
        if (warp.warp.top.it > 0) and (warp.warp.top.it%self._Plot == 0):
            return self.GetRenderScripts().keys()
        else:
            return []

    #------------------------------------------------------------------------
    def Advance(self):
        """
        Advance the simulation n steps, where n can be passed on the
        command line via --step option or by calling SetStepInterval.
        """
        warp.step(self._Step)
        return

    #------------------------------------------------------------------------
    def Continue(self):
        """
        Return False when the run is finished, and True otherwise.
        This implementation will run until the n'th simulation step
        where n can be passed on the command line via the --stop
        option or by calling SetStopIteration.
        """
        return warp.warp.top.it < self._Stop


    #------------------------------------------------------------------------
    def Finalize(self):
        """
        This is called as part of normal shutdown. Clean up code can be
        insrted here. This implementation does nothing.
        """
        return

    #------------------------------------------------------------------------
    def Initialize(self):
        """
        This will be called as part of start up. This method should be
        overriden to configure the simulation initial condition etc.

        This implementation takes n steps, where n may be passed on the
        command line via the --start option, or by calling SetStartIteration.
        """
        warp.step(self._Start)
        return

    # the following helper functions support the default implementation
    # but are not part of the interface that a user needs to implement
    #-------------------------------------------------------------------------
    def AddRenderScript(self, key, fileName):
        """
        Append a render script to the dictionary. This should be called
        any number of times before Initialize runs.
        """
        if self._RenderScripts is None:
            self._RenderScripts = {}
        self._RenderScripts[key] = fileName
        return

    #-------------------------------------------------------------------------
    def SetStartIteration(self, it):
        """
        Set the number of simulation steps to take before visualization
        """
        self._StartIteration = it
        return

    #-------------------------------------------------------------------------
    def SetStopIteration(self, it):
        """
        Set the simulation iteration to stop at.
        """
        self._StopIteration = it
        return

    #-------------------------------------------------------------------------
    def SetStepInterval(self, it):
        """
        Set the number of iterations the simulation advances per update.
        """
        self._Step = it
        return

    #-------------------------------------------------------------------------
    def SetPlotInterval(self, it):
        """
        Set the number of iterations between visualizations.
        """
        self._Plot = it
        return
