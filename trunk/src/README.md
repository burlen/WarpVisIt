#WarpVisIt in-situ visualization scripts.
This porject contains code to use VisIt in-situ with the Warp simulation.

##Project Organization
* scripts
    * This directory contains simulation and visualization scripts
* site
    * This directory contains site specific configuration.

##Install
TODO -- python setup.py install?? CMake?? Some other??

##Running
The user's shell enviornment must be configured with both Warp and VisIt in the PATH. This is system specific, eg. on NERSC systems there will be modules for this. The shell command WarpVisItEnv.sh can be used to configure the VisIt enviornment.

Once the environment is configured runs are started using the WarpVisIt.sh command. The user has to at least provide a Warp configuration script, and may optionally provide a VisIt script for non-interactive runs. For more information see the output of WarpVisIt.sh --help.

###Warp Script
The user must provide a script to configure the simulation with the following functions defined:

* void Initialize()
    * This function contains the Warp specific code to configure the simulation, and initialize MPI and Warp.
* void Finalize()
    * This function contains code to shut the simulation down and finalize MPI and Warp.
* bool Finished()
    * This function returns true when the simulation is finished, ie taken the desired number of steps.

See the scripts directory for examples.

###VisIt Script
When running interactively the user can instead open the Sim2 file in the VisIt GUI and configure visualizations interactively. When running non-interactively, the user must provide a script that tells VisIt what to do with the following functions defined:

* void Render(string sim2FileName)
    * This function is passed the sim2 file and contains the code to drive VisIt rendering.

See the scripts directory for examples.

##Developer documentation
###Sources
* WarpVisItEnv.sh
    * A helper script that given the path to a VisIt install, will setup the shell environment.

* WarpVisIt.sh
    * The script used to launch runs.

* WarpVisItMain.py
    * The program's main routine.

* WarpVisItMPIMain.sh
    * A wrapper around the program's main, it is needed to pre-load OpenGL symbols.

* WarpVisItSimV2Db.py
    * VisIt simV2 database code callbacks.

* WarpVisItEngine.py
    * code that drives the VisIt engine.

* WarpVisItViewer.py
    * code that drives the VisIt viewer. This comes into play only for non-interactive runs.

* WarpVisItUtil.py
    * a collection of helper routines.

###Environment variables
The user supplied scripts are passed into the runtime by the following environment variables:

* WARP_SCRIPT
    * This variable points to a script contains the Warp simulation setup code.
* VISIT_SCRIPT
    * If provided this variable points to a sript that tells VisIt what to do. If this variable is provided then the run will be non-interactive. If this variable is not provided then the run will be interactive and the user must connect to the simulation from the VisIt GUI. This is accomplished by opening the sim2 file.
* SIM2_FILE
    * If provided the variable points to where VisIt should place the .sim2 file that allows it's GUI to connect to the simulation interactively. If not provided the file will be called Warp.sim2 and placed in the current working directory(CWD).
