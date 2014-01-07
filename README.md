#WarpVisIt in-situ visualization scripts.
This project contains code to use VisIt in-situ with the Warp simulation.

##Install
TODO -- python setup.py install?? CMake?? Some other??

##Running
The *WarpVisIt.sh* is used to launch runs. The user has to at least provide a simulation script, and may optionally provide any number of rendering scripts for non-interactive runs. For more information see the output of WarpVisIt.sh --help.

###Interactive Example
./WarpVisIt.sh -np=8  --warp-script=/work/warp-visit/WarpVisIt/scripts/cigar/cigar.py --visit-install=/work/warp-visit/visit/visit-install/ --sim2-file=/work/warp-visit/visit/Warp.sim2 --script-dir=/work/warp-visit/WarpVisIt/scripts/cigar --interactive

###Non-interactive Example
./WarpVisIt.sh -np=8  --warp-script=/work/warp-visit/WarpVisIt/scripts/cigar/cigar.py --visit-install=/work/warp-visit/visit/visit-install/ --sim2-file=/work/warp-visit/visit/Warp.sim2 --script-dir=/work/warp-visit/WarpVisIt/scripts/cigar

###Simulation script
The user must provide a script to configure the simulation and visualizaitons. The API is as follows:

* *LoadRenderScripts(scriptRoot)* -- Function to load rendering scripts and assign them unique names.
* *GetActiveRenderScripts()* -- Function that returns a list of render script names to render each iteration.
* *Advance()* -- Function that when called advances the simulation by one step.
* *Continue()* -- Function that returns True while the simulation should scontinue to advance.
* *Finalize()* -- Function to cleaup and release resources as the simulation ends.
* *Initialize()* -- Funciton to alocate resources, configure initial conditions etc, as the simulation starts.

 Any number of additional visualization scripts may be provided, they are loaded by the LoadRenderScripts function above. See the scripts directory for examples.

##Developer documentation
###Project organization
* *scripts* -- This directory contains simulation and visualization scripts
* *site* -- This directory contains site specific configuration.

###Class documentation
* *WarpVisItEnv.sh* -- A helper script that given the path to a VisIt install, will setup the shell environment.
* *WarpVisIt.sh* -- The script used to launch runs.
* *WarpVisItMain.py* -- The program's main routine.
* *WarpVisItMPIMain.sh* -- A wrapper around the program's main, it is needed to pre-load OpenGL symbols.
* *WarpVisItSimV2Db.py* -- VisIt simV2 database code callbacks.
* *WarpVisItEngine.py* -- code that drives the VisIt engine.
* *WarpVisItCLI.py* -- code that drives the CLI. This comes into play only for non-interactive runs.
* *WarpVisItUtil.py* -- a collection of helper routines.

###Environment variables
The user's shell enviornment must be configured with both Warp and VisIt in the PATH. This is system specific, eg. on NERSC systems there will be modules for this. The shell command WarpVisItEnv.sh can be used to configure the VisIt enviornment.
The user supplied scripts are passed into the runtime by the following environment variables:

* *WARPVISIT_WARP_SCRIPT* -- This variable points to a script contains the Warp simulation setup code.
* *ARPVISIT_SCRIPT_DIR* -- This points to the directory containing any user provided scripts.
* *WARPVISIT_SIM2_FILE* -- If provided the variable points to where VisIt should place the .sim2 file that allows it's GUI to connect to the simulation interactively. It defaults to the current working directory(CWD)
