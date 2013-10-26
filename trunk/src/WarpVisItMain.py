#!python
import os
import imp
from WarpVisItUtil import VisItEnv
env = VisItEnv()
from WarpVisItUtil import pError
from WarpVisItEngine import WarpVisItEngine
from WarpVisItCLI import WarpVisItCLI

def main(argv=None):

    # load user specified simulation configuration.
    # this script is required.
    warpScriptFileName = os.getenv('WARP_SCRIPT')
    if not warpScriptFileName:
        pError('WARP_SCRIPT environment variable was not set.')
        return -1
    if not os.path.isfile(warpScriptFileName):
        pError('The file WARP_SCRIPT points to (%s) was not found'%simConfigFileName)
        return -1
    warpScript = imp.load_source('warpScript', warpScriptFileName)
    warpScript.Initialize()

    # Initialize the VisIt engine. The engine runs in the
    # simulation address space.
    engine = WarpVisItEngine()

    simFileName = os.getenv('SIM2_FILE')
    if simFileName:
        engine.SetSimFile(simFileName)

    engine.Initalize()

    # load user specified rendering configuration (optional)
    # if it's provided we'll fork a CLI and Viewer. The CLI
    # will then take control of the run otherwise CLI does
    # nothing and the Viewer is assumed external.
    cli = WarpVisItCLI()
    cli.SetSimFile(engine.GetSimFile())

    visitScriptFileName = os.getenv('VISIT_SCRIPT')
    if visitScriptFileName:
        if not os.path.isfile(visitScriptFileName):
            pError('VISIT_SCRIPT %s not found'%visitScriptFileName)
            return -1
        visitScript = imp.load_source('visitScript', visitScriptFileName)

        cli.SetRenderCallback(visitScript.Render)
        cli.SetShouldRenderCallback(warpScript.UpdatePlots)

    cli.Initialize()

    # the engine now runs until either the desired number of
    # steps is reached, or the cli tells him to stop.
    while ((not warpScript.Finished())
            and engine.EventLoop() and cli.EventLoop()):
        pass

    # shut the engine and cli down
    engine.Finalize()
    cli.Finalize()

    return 0

if __name__ == "__main__":
    main()
