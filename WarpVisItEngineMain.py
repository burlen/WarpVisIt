#!/bin/python
import os
import sys
# hide command line from warp
argv = sys.argv
sys.argv = [argv[0]]
#
import argparse
import time
from WarpVisItUtil import VisItEnv
env = VisItEnv()
from WarpVisItUtil import pError,pDebug,pStatus
from WarpVisItEngine import WarpVisItEngine
from WarpVisItSimulation import WarpVisItSimulation
from WarpVisItSimulationFactory import WarpVisItSimulationFactory

#-----------------------------------------------------------------------------
if __name__ == "__main__":
    pStatus('WarpVisItEngineMain started')

    # parse command line args
    ap = argparse.ArgumentParser(usage=argparse.SUPPRESS,prog='',add_help=False)
    ap.add_argument('--help',default=False,action='store_true')
    opts = vars(ap.parse_known_args(argv)[0])
    if opts['help']:
        pStatus('\nWarpVisItEngineMain\nUsage:\n\nWarpVisItEngine\n%s\nWarpVisItSimulation\n%s\nWarpVisItSimulationFactory\n%s'%(
            WarpVisItEngine.GetCommandLineHelp(),
            WarpVisItSimulation.GetCommandLineHelp(),
            WarpVisItSimulationFactory.GetCommandLineHelp()))
        sys.exit(0)

    # use the factory to create user's simualtion
    factory = WarpVisItSimulationFactory(argv)
    simulation = factory.CreateSimulation()
    simulation.Initialize()

    # create a visit engine
    engine = WarpVisItEngine(argv)
    engine.Initalize()
    engine.SetSimulation(simulation)

    # run
    status = engine.EventLoop()

    # shut down engine and simulation
    engine.Finalize()
    simulation.Finalize()

    pStatus('WarpVisItEngineMain finished')
    sys.exit(status)
