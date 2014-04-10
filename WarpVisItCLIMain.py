#!python
import os
import sys
import argparse
import time
from WarpVisItUtil import VisItEnv
env = VisItEnv()
from WarpVisItUtil import pError
from WarpVisItUtil import pDebug
from WarpVisItCLI import WarpVisItCLI

#-----------------------------------------------------------------------------
if __name__ == "__main__":

    # parse command line args
    ap = argparse.ArgumentParser(usage=argparse.SUPPRESS,prog='',add_help=False)
    ap.add_argument('--help',default=False,action='store_true')
    opts = vars(ap.parse_known_args(sys.argv)[0])
    if opts['help']:
        sys.stderr.write('WarpVisItCLIMain\nUsage:\n\n%s\n\n'%(WarpVisItCLI.GetCommandLineHelp()))
        sys.exit(0)

    # create and run CLI/viewer
    cli = WarpVisItCLI(sys.argv)
    cli.Initialize()
    status = cli.EventLoop()
    cli.Finalize()

    sys.exit(status)
