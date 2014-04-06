#!python
import os
import sys
import imp
import getopt
import time
from WarpVisItUtil import VisItEnv
env = VisItEnv()
from WarpVisItUtil import pError
from WarpVisItUtil import pDebug
#from WarpVisItEngine import WarpVisItEngine
#from WarpVisItCLI import WarpVisItCLI, CLIMain

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
    from WarpVisItEngine import WarpVisItEngine
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
    #engine.Initalize(viewerOpts=['-nowin','-debug','1'],engineOpts='-debug 1')

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
    simFile = os.path.abspath(os.getenv('WARPVISIT_SIM2_FILE'))
    scriptRoot = os.getenv('WARPVISIT_SCRIPT_DIR')
    interactive = bool(os.getenv('WARPVISIT_INTERACTIVE'))
    cli = bool(os.getenv('WARPVISIT_CLI'))

    # then command line arguments
    opts, args = getopt.getopt(sys.argv,'',
        ['help','warp-script=','script-dir=','sim-file=','interactive','cli'])

    for opt, arg in opts:
        if opt == '--warp-script':
            warpScript = arg
        elif opt == '--script-dir':
            scriptRoot = arg
        elif opt == '--sim-file':
            simFile = os.path.abspath(arg)
        elif opt == '--interactive':
            interactive = bool(arg)
        elif opt == '--cli':
            cli = bool(arg)
        else:
            sys.stderr.write(
                """
                command line arguments and environment variable:
                    --warp-script : WARPVISIT_WARP_SCRIPT : required
                    --script-dir  : WARPVISIT_SCRIPT_DIR  : optional
                    --sim-file    : WARPVISIT_SIM2_FILE   : optional
                    --interactive : WARPVISIT_INTERACTIVE : optional
                    --cli         : WARPVISIT_CLI         : optional
                """)
            sys.exit(0)

    # then validate
    if not cli and not warpScript:
        pError('WARP_SCRIPT environment variable was not set.')
        sys.exit(-1)

    if not cli and not os.path.isfile(warpScript):
        pError('The file WARP_SCRIPT points to (%s) was not found'%warpScript)
        sys.exit(-1)

    if scriptRoot and not os.path.isdir(scriptRoot):
        pError('The path SCRIPT_DIR points to (%s) was not found'%scriptRoot)
        sys.exit(-1)

    if not cli and not simFile:
        pError('running a CLI and no SIM2_FILE was given.')
        sys.exit(-1)

    # finally start
    if cli:
        # CLI/viewer
        from WarpVisItCLI import WarpVisItCLI, CLIMain
        timeout = 1000
        simFileFound=False
        n = timeout/10
        i = 0
        while i<n:
            if os.path.isfile(simFile):
                simFileFound=True
                sys.stderr.write('found %s'%(simFile))
                break
            else:
                sys.stderr.write('.')
                time.sleep(10)
            i += 1
        sys.stderr.write('\n')

        if simFileFound:
            status = CLIMain(simFile)
        else:
            sys.stderr.write('Error: failed to find %s'%(simFile))
            status = -1

        sys.exit(status)

    else:
        # engine/simulation
        status = main(warpScript, simFile, scriptRoot, interactive)
        sys.exit(status)
