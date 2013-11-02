#!python
import os
import sys
import imp
import getopt
from WarpVisItUtil import VisItEnv
env = VisItEnv()
from WarpVisItUtil import pError
from WarpVisItUtil import pDebug
from WarpVisItEngine import WarpVisItEngine
from WarpVisItCLI import WarpVisItCLI

#-----------------------------------------------------------------------------
def main(warpScriptFileName='', simFileName='', scriptRoot='',
    interactive=True):
    """
    The program main routine. It takes the following arguments:

        warpScript -- a path/file that has user provided code (required)
        sim2FileName -- a path/file indicating where to write sim2 file (optional)
        scriptRoot -- a path to where scripts can be found (optional)
        interactive -- a flag indicating interactive or headless run (optional)
    """
    pDebug('main %s %s %s %s'%(
        warpScriptFileName, simFileName, scriptRoot, str(interactive)))

    # load user specified simulation configuration.
    # this script is required.
    warpScript = imp.load_source('warpScript', warpScriptFileName)
    warpScript.Initialize()

    # Initialize the VisIt engine. The engine runs in the
    # simulation address space.
    engine = WarpVisItEngine()
    engine.SetInteractive(interactive)
    engine.SetSimFile(simFileName)
    engine.SetStepCallback(warpScript.Advance)
    engine.SetContinueCallback(warpScript.Continue)
    engine.SetActiveRenderScriptsCallback(warpScript.GetActiveRenderScripts)
    engine.SetRenderScripts(warpScript.LoadRenderScripts(scriptRoot))
    engine.Initalize()

    # the engine now runs until either the desired number of
    # steps is reached, or the cli tells him to stop.
    engine.EventLoop()

    # shut the engine and cli down
    engine.Finalize()
    warpScript.Finalize()

    return 0

#-----------------------------------------------------------------------------
if __name__ == "__main__":

    # environment variables first
    warpScript = os.getenv('WARPVISIT_WARP_SCRIPT')
    simFile = os.getenv('WARPVISIT_SIM2_FILE')
    scriptRoot = os.getenv('WARPVISIT_SCRIPT_DIR')
    interactive = bool(os.getenv('WARPVISIT_INTERACTIVE'))

    # then command line arguments
    opts, args = getopt.getopt(sys.argv,'',
        ['help','warp-script=','script-dir=','sim-file=','interactive='])

    for opt, arg in opts:
        if opt == '--warp-script':
            warpScript = arg
        elif opt == '--script-dir':
            scriptRoot = arg
        elif opt == '--sim-file':
            simFile = arg
        elif interactive == '--interactive':
            interactive = bool(arg)
        else:
            sys.stderr.write(
                """
                command line arguments and environment variable:
                    --warp-script : WARPVISIT_WARP_SCRIPT : required
                    --script-dir  : WARPVISIT_SCRIPT_DIR  : optional
                    --sim-file    : WARPVISIT_SIM2_FILE   : optional
                    --interactive : WARPVISIT_INTERACTIVE : optional
                """)
            sys.exit(0)

    # then validate
    if not warpScript:
        pError('WARP_SCRIPT environment variable was not set.')
        sys.exit(-1)
    if not os.path.isfile(warpScript):
        pError('The file WARP_SCRIPT points to (%s) was not found'%warpScript)
        sys.exit(-1)

    if scriptRoot and not os.path.isdir(scriptRoot):
        pError('The path SCRIPT_DIR points to (%s) was not found'%scriptRoot)
        sys.exit(-1)

    # finally run with it
    status = main(warpScript, simFile, scriptRoot, interactive)
    sys.exit(status)

#    # load user specified rendering configuration (optional)
#    # if it's provided we'll fork a CLI and Viewer. The CLI
#    # will then take control of the run otherwise CLI does
#    # nothing and the Viewer is assumed external.
#    cli = WarpVisItCLI()
#    cli.SetSimFile(engine.GetSimFile())
#
#    visitScriptFileName = os.getenv('VISIT_SCRIPT')
#    if visitScriptFileName:
#        if not os.path.isfile(visitScriptFileName):
#            pError('VISIT_SCRIPT %s not found'%visitScriptFileName)
#            return -1
#        visitScript = imp.load_source('visitScript', visitScriptFileName)
#
#        cli.SetRenderCallback(visitScript.Render)
#        cli.SetShouldRenderCallback(warpScript.UpdatePlots)
#
#    cli.Initialize()

