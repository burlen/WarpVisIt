#!/bin/python
import os
import sys
# hide command line from warp
import warpoptions
warpoptions.ignoreUnknownArgs = True
warpoptions.quietImport = True
warpoptions.init_parser()

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
    opts = vars(ap.parse_known_args(sys.argv)[0])
    if opts['help']:
        pStatus('\nWarpVisItEngineMain\nUsage:\n\nWarp\n%s\nWarpVisItEngine\n%s\nWarpVisItSimulation\n%s\nWarpVisItSimulationFactory\n%s'%(
            warpoptions.warpoptionsstr(),
            WarpVisItEngine.GetCommandLineHelp(),
            WarpVisItSimulation.GetCommandLineHelp(),
            WarpVisItSimulationFactory.GetCommandLineHelp()))
        sys.exit(0)

    # use the factory to create user's simualtion
    factory = WarpVisItSimulationFactory(sys.argv)
    simulation = factory.CreateSimulation()
    simulation.Initialize()

    # create a visit engine
    engine = WarpVisItEngine(sys.argv)
    engine.Initalize()
    engine.SetSimulation(simulation)

    # run
    status = engine.EventLoop()

    # shut down engine and simulation
    engine.Finalize()
    simulation.Finalize()

    pStatus('WarpVisItEngineMain finished')
    sys.exit(status)
