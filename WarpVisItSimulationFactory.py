from WarpVisItUtil import pError
from WarpVisItUtil import pDebug
import os
import argparse
import sys
import imp

class WarpVisItSimulationFactory:
    """
    Create an instance of WarpVisItSimulation using the
    user supplied factory function. This function is expected
    to be at file scope of the given script and named
    NewWarpVisItSimulation
    """
    #-------------------------------------------------------------------------
    def __init__(self, args=[]):
        """  """
        self.__Args = args
        self.__FactoryScript = os.getenv('WARPVISIT_FACTORY_SCRIPT')
        if not self.__FactoryScript:
            self.__FactoryScript = ''

        # parse command line args
        ap = argparse.ArgumentParser(usage=argparse.SUPPRESS,prog='WarpVisItSimulationFactory',add_help=False)
        ap.add_argument('--factory-script',type=str,default=self.__FactoryScript)
        ap.add_argument('--script-dir',type=str,default=None)
        opts = vars(ap.parse_known_args(args)[0])
        self.__FactoryScript = os.path.abspath(opts['factory_script'])

        return

    #-------------------------------------------------------------------------
    def SetFactoryScript(self, scriptFile):
        """
        Set the path to the script that contains the
        file scope factory function named NewWarpVisItSimulation
        """
        self.__FactoryScript = scriptFile

    #-------------------------------------------------------------------------
    @staticmethod
    def GetCommandLineHelp():
        """
        Return list of command line options and environment vars
        """
        return
        "--factory-script : WARPVISIT_FACTORY_SCRIPT : path to WarpVisItSimulation factory\n"

    #-------------------------------------------------------------------------
    def CreateSimulation(self):
        """
        Create a WarpVisItSimulation instance using the user
        supplied factory function.
        """
        if not os.path.isfile(self.__FactoryScript):
            pError('Failed to locate factory script (%s)'%self.__FactoryScript)
            raise RuntimeError('Invalid script file')

        # load the script containing the factory function
        factory = imp.load_source('warpScript', self.__FactoryScript)
        # use the factory to create an interface we can use to
        # which ever simulation the user is running
        simulation = factory.NewWarpVisItSimulation(self.__Args)
        return simulation
