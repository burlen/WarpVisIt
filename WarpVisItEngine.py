import os
import sys
import argparse
import simV2
import WarpVisItSimV2Db
from WarpVisItUtil import pError, pDebug, pStatus, VisItEnv, getEnvVar
import traceback
import parallel


#-----------------------------------------------------------------------------
class WarpVisItEngine:
    """
    The WarpVisItEngine manages engine state and exposes
    engine API specifically to the Warp simulation.
    """

    #-------------------------------------------------------------------------
    def __init__(self, args):
        """ """
        self._CommRank = parallel.get_rank()
        self._CommSize = parallel.number_of_PE()
        self._Env = VisItEnv()
        self._SimFile = getEnvVar('WARPVISIT_SIM_FILE',str,'WarpVisIt.sim2')
        self._TraceFile = getEnvVar('WARPVISIT_TRACE_FILE',bool,False)
        self._Simulation = None
        self._ShutdownRequested = False
        self._Connected = False
        self._Synchronous = 0
        self._Behavior = BatchBehavior(self)
        self._UpdateVisItGUI = True
        self._VisItControlStepping = True
        self._VisItBlockingComm = True
        self._CommandQueue = []
        self._EngineOpts = None
        self._StepCount = 0
        self._ProbeMem = 0

        # parse command line args
        interact = getEnvVar('WARPVISIT_MODE_INTERACT',bool,False)
        monitor = getEnvVar('WARPVISIT_MODE_MONITOR',bool,False)
        batch = getEnvVar('WARPVISIT_MODE_BATCH',bool,False)

        ap = argparse.ArgumentParser(usage=argparse.SUPPRESS,prog='WarpVisItEngine',add_help=False)
        ap.add_argument('--sim-file',type=str,default=self._SimFile)
        ap.add_argument('--trace-file',default=self._TraceFile, action='store_true')
        ap.add_argument('--interact',default=interact, action='store_true')
        ap.add_argument('--monitor',default=monitor, action='store_true')
        ap.add_argument('--prompt',default=False, action='store_true')
        ap.add_argument('--batch',default=batch, action='store_true')
        ap.add_argument('--engine-opts',type=str,default=None)
        ap.add_argument('--probe-mem',type=int,default=self._ProbeMem)
        opts = vars(ap.parse_known_args(args)[0])
        self._SimFile = os.path.abspath(opts['sim_file'])
        self._TraceFile = opts['trace_file']
        interact = opts['interact']
        monitor = opts['monitor']
        prompt = opts['prompt']
        batch = opts['batch']
        self._EngineOpts = opts['engine_opts']
        self._ProbeMem = opts['probe_mem']

        if sum((interact,monitor,prompt,batch))>1:
            raise RuntimeError('--interact, --monitor, --batch are mutually exclusive')
        if interact:
            self._Behavior = InteractBehavior(self)
        elif monitor:
            self._Behavior = MonitorBehavior(self)
        elif prompt:
            self._Behavior = PromptBehavior(self)
        else:
            self._Behavior = BatchBehavior(self)

        # initialize libsim
        simV2.VisItSetDirectory(self._Env.GetRoot())
        simV2.VisItSetupEnvironment()
        simV2.VisItSetParallelRank(self._CommRank)
        simV2.VisItSetParallel(self._CommSize > 1)
        simV2.VisItSetBroadcastIntFunction(broadcastInt)
        simV2.VisItSetBroadcastStringFunction(broadcastString)

        return


    #-------------------------------------------------------------------------
    @staticmethod
    def GetCommandLineHelp():
        """
        Return list of command line options and environment vars
        """
        return ("--sim-file : WARPVISIT_SIM_FILE : Path to write sim2 file\n"
           "--trace-file : WARPVISIT_TRACE_FILE : Path to write trace file\n"
           "--interact : WARPVISIT_MODE_INTERACT : Wait for VisIt and run with VisIt control\n"
           "--monitor : WARPVISIT_MODE_MONITOR : Listen for while running, continue without VisIt\n"
           "--batch : WARPVISIT_MODE_BATCH : Run without VisIt\n"
           "--engine-opts : : Command line options to pass into the VisIt engine\n"
           "--probe-mem : : Number of sim steps between memory probes\n")


    #-------------------------------------------------------------------------
    def SetSimulation(self, simulation):
        """Set the simulation object"""
        self._Simulation = simulation

    #-------------------------------------------------------------------------
    def RequestShutdown(self):
        """Request shutdown"""
        self._VisItControlStepping = False
        self._VisItBlockingComm = False
        self._UpdateVisItGUI = False
        self._ShutdownRequested = True

    #-------------------------------------------------------------------------
    def SetVisItControlStepping(self, v):
        """Select between interactive/non-interactive simulation stepping"""
        self._VisItControlStepping = bool(v)
        self._VisItBlockingComm = bool(v)

    #-------------------------------------------------------------------------
    def GetVisItControlStepping(self, v):
        """Select between interactive/non-interactive simulation stepping"""
        return self._VisItControlStepping

    #-------------------------------------------------------------------------
    def SetSynchronous(self, v):
        """
        Tracks wether we are currently rendering. we need to prevent simulation
        advances while rendering, else we may mix data from multiple time steps
        in the same plots.
        """
        if bool(v):
            self._Synchronous += 1
        else:
            self._Synchronous -= 1
        if self._Synchronous < 0:
            self._Synchronous = 0

    #-------------------------------------------------------------------------
    def GetSynchronous(self):
        """True when VisIt is actively rendering"""
        return self._Synchronous > 0

    #-------------------------------------------------------------------------
    def SetSimFile(self, fileName):
        """
        Set the .sim2 file name. If one is not set the default is
        '/cwd/Warp.sim2'
        """
        self._SimFile = fileName

    #-------------------------------------------------------------------------
    def GetSimFile(self):
        """Return the .sim2 file name"""
        return self._SimFile

    #-------------------------------------------------------------------------
    def SetInteractive(self, interactive):
        """Set true if the simulation should connect to VisIt"""
        self._Interactive = bool(interactive)
        self._UpdateVisItGUI = bool(interactive)

    #-------------------------------------------------------------------------
    def Initalize(self, engineOpts=''):
        """
        Perform the initial setup of VisIt to include the libsim module
        with the simulation.
        """
        if __debug__: pDebug('WarpVisItEngine::Initialize')

        simV2.VisItSetDirectory(self._Env.GetRoot())

        if self._TraceFile:
            simV2.VisItOpenTraceFile('%d.trace'%(self._CommRank))

        if self._EngineOpts:
            if __debug__: pDebug('Using options %s'%(self._EngineOpts))
            simV2.VisItSetOptions(self._EngineOpts);

        if self._CommRank == 0:
            if not simV2.VisItInitializeSocketAndDumpSimFile(
                  'WarpVisIt', None, os.getcwd(), None, None, self._SimFile):
                pError('VisIt initialization failed')
                return False

        self._Behavior.Initialize()

        return True

    #-------------------------------------------------------------------------
    def Finalize(self):
        """
        Finalize VisIt libsim, cleanup, etc
        """
        if self._TraceFile:
            simV2.VisItCloseTraceFile()

        # rm this file as if its left around
        # it will cause next cli to fail
        if self._CommRank == 0:
            if os.path.isfile(self._SimFile):
                os.unlink(self._SimFile)

        return

    #-------------------------------------------------------------------------
    def EventLoop(self):
        """
        Run the VisIt controlled event loop. Return False on error
        """
        doError = lambda e: e<0
        doUpdate = lambda e: e==0
        doConnect = lambda e: e==1
        doInternal = lambda e: e==2
        masterRank = lambda r: r==0
        onErrorContinue = lambda : self._Behavior=='Monitor'

        if __debug__: pDebug('WarpVisItEngine::EventLoop')

        while 1:
            # test or wait for incomming communication.
            if masterRank(self._CommRank):
                event = simV2.VisItDetectInput(self._VisItBlockingComm, -1)
                parallel.broadcast(event)
            else:
                event = parallel.broadcast(0)

            # process incomming communication
            if __debug__: pDebug('event=%d'%(event))

            # internal comm error
            if doError(event):
                # the VisIt viewer disconnected.
                pError('An communication error was detected')
                self._Behavior.Initialize()
                if onErrorContinue():
                    continue
                else:
                    return False

            # simulate
            elif doUpdate(event):
                # this is where we can safely do things.
                # it occurs only when VisIt is not in the
                # middle of processing its own commands.
                if __debug__: pDebug('update')
                self._Behavior.Update()

            # internal connection
            elif doConnect(event):
                if simV2.VisItAttemptToCompleteConnection() == 1:
                    pStatus('WarpVisItEngine connected')
                    self.ConnectLibsim()
                    self._Behavior.Connect()
                    self._Behavior.Update()
                else:
                    pError('Connection failed')
                    self._Behavior.Initialize()
                    if onErrorContinue():
                        continue
                    else:
                        return False

            # internal comm
            elif doInternal(event):
                if __debug__: pDebug('internal')
                if not self.CommunicateLibsim():
                    pError('Error while processing internal commands')
                    self._Behavior.Initialize()
                    if onErrorContinue():
                        continue
                    else:
                        return False

            # bad value
            else:
                pError('unknown command %d'%(event))

            # check for asynchronous shutdown request
            if self._ShutdownRequested:
                pStatus('WarpVisItEngine exiting the event loop')
                simV2.VisItDisconnect()
                self._Connected = False
                return True

            continue

        return False

    #------------------------------------------------------------------------
    def ConnectLibsim(self):
        """
        Install database callbacks in libsim
        """
        if __debug__: pDebug('WarpVisItEngine::InitializeLibsim')
        simV2.VisItSetGetMetaData(WarpVisItSimV2Db.getMetaData, self)
        simV2.VisItSetGetMesh(WarpVisItSimV2Db.getMesh, self)
        simV2.VisItSetGetVariable(WarpVisItSimV2Db.getVar, self)
        simV2.VisItSetGetDomainList(WarpVisItSimV2Db.getDomains, self)
        simV2.VisItSetCommandCallback(commandCallback, self)
        simV2.VisItSetSlaveProcessCallback(slaveProcessCallback)
        self._Connected = True
        self._VisItBlockingComm = True
        return

    #-------------------------------------------------------------------------
    def CommunicateLibsim(self):
        """ProcessCommand internal VisIt commands."""
        masterRank = lambda r: r==0
        COMMAND_PROCESS = 0
        COMMAND_SUCCESS = 1
        COMMAND_FAILURE = 2

        if __debug__: pDebug('WarpVisItEngine::ProcessCommands')

        if masterRank(self._CommRank):
            if simV2.VisItProcessEngineCommand():
                if __debug__: pDebug('master success')
                broadcastSlaveCommand(COMMAND_SUCCESS)
                return True

            else:
                if __debug__: pDebug('master discconnect')
                broadcastSlaveCommand(COMMAND_FAILURE)
                return False

        else:
            while True:
                command = broadcastSlaveCommand()
                if __debug__: pDebug('slave command %i'%(command))

                if command == COMMAND_PROCESS:
                    if __debug__: pDebug('slave process command')
                    simV2.VisItProcessEngineCommand()

                elif command == COMMAND_SUCCESS:
                    if __debug__: pDebug('slave success')
                    return True

                elif command == COMMAND_FAILURE:
                    if __debug__: pDebug('slave disconnect')
                    return False

                else:
                    pError('Bad command %i'%(command))

        return False

    #-------------------------------------------------------------------------
    def ProcessCommand(self, cmd):
        """
        process simulation control commands.

        endSyn    : internal use only, for synchronziation between sim and vis

        """
        if __debug__: pDebug('WarpVisItEngine::ProcessCommand %s'%(cmd))

        self._CommandQueue.append(cmd)

        # queue commands while in synchronous mode
        if self.GetSynchronous():
            if cmd == 'endSyn':
                self.SetSynchronous(False)
                self._CommandQueue.pop()
            else:
                if __debug__: pDebug('defered %s'%(cmd))
                return

        # now in asynchronous mode
        # process queued commands
        ncmds = len(self._CommandQueue)
        i = 0
        while (i < ncmds):
            qcmd = self._CommandQueue.pop()
            self._Behavior.Process(qcmd)
            i += 1
        return

    #-------------------------------------------------------------------------
    def PushRenderScripts(self):
        """Render simulation data"""
        if __debug__: pDebug('WarpVisItEngine::Render')

        activeScripts = self._Simulation.GetActiveRenderScripts()
        if len(activeScripts):
            scripts = self._Simulation.GetRenderScripts()
            # push rendering scripts for execution
            # while rendering sim must not advance
            # enable synchronous mode for each script
            # and append a command to disable it to
            # the end of each.
            for script in activeScripts:
                if __debug__: pDebug('Rendering %s'%(script))
                self.SetSynchronous(True)
                source = scripts[script]
                source += ("\nfrom visit import visit\n"
                           "visit.SendSimulationCommand('localhost', '%s', 'endSyn')\n")%(
                          self._SimFile)
                simV2.VisItExecuteCommand(source)

        return True

    #-------------------------------------------------------------------------
    def StepSimulation(self):
        """Drive the simulation through one or more steps"""
        if __debug__: pDebug('WarpVisItEngine::StepSimulation')

        self._Simulation.Advance()
        simV2.VisItTimeStepChanged()
        self._StepCount += 1
        return True

    #-------------------------------------------------------------------------
    def ProbeMemory(self,force=False):
        """ """
        # report memory consumption
        if force or (self._ProbeMem > 0) and ((self._StepCount%self._ProbeMem) == 0):
            ret,vm,rss = simV2.VisItGetMemory()
            vmn = parallel.globalmin(vm)
            vmx = parallel.globalmax(vm)
            vme = parallel.globalave(vm)
            pStatus('MemUse=%g %g %g'%(vmn,vme,vmx))
        return

    #-------------------------------------------------------------------------
    def Disconnect(self):
        """
        Handle VisIt disconnect
        """
        if self._Connected:
            simV2.VisItDisconnect()
            self._Connected = False
        return









##############################################################################
class InteractBehavior(object):
    """Implements interactive behavior"""
    #-------------------------------------------------------------------------
    def __init__(self, engine):
        """ """
        self.__Engine = engine
        return

    #-------------------------------------------------------------------------
    def Initialize(self):
        """
        Initialize the object for interact mode
        """
        self.__Engine.Disconnect()
        self.__Engine._UpdateVisItGUI = False
        self.__Engine._VisItControlStepping = True
        self.__Engine._VisItBlockingComm = True
        self.__Engine._Synchronous = 0
        self.__Engine._CommandQueue = []
        return

    #-------------------------------------------------------------------------
    def Connect(self):
        """
        Connect the object for interact mode
        """
        self.__Engine._UpdateVisItGUI = True
        self.__Engine._VisItControlStepping = True
        self.__Engine._VisItBlockingComm = True
        self.__Engine._Synchronous = 0
        self.__Engine._CommandQueue = []
        return

    #-------------------------------------------------------------------------
    def Update(self):
        """
        Respond to request for render in interact mode
        """
        # interactive update
        if self.__Engine._Connected:
            # render plot for this sim time step
            if self.__Engine._UpdateVisItGUI:
                # trigger render of plot defined in the GUI
                simV2.VisItUpdatePlots()

            # advance the simulation if auto stepping
            if not self.__Engine._VisItControlStepping:
                self.__Engine.StepSimulation()
        return

    #-------------------------------------------------------------------------
    def Process(self, qcmd):
        """
        process interact control commands.
        """
        if (qcmd == 'end'):
            self.__Engine.RequestShutdown()

        elif qcmd == 'pause':
            pStatus('pause')
            self.__Engine._VisItBlockingComm = True
            self.__Engine._VisItControlStepping = True
            self.__Engine._UpdateVisItGUI = True
            self.Update()

        elif qcmd == 'step':
            pStatus('step')
            self.__Engine._VisItBlockingComm = True
            self.__Engine._VisItControlStepping = True
            self.__Engine._UpdateVisItGUI = True
            self.__Engine.StepSimulation()
            self.Update()

        elif qcmd == 'run':
            pStatus('run')
            self.__Engine._VisItBlockingComm = False
            self.__Engine._VisItControlStepping = False
            self.__Engine._UpdateVisItGUI = True

        elif qcmd == 'continue':
            pStatus('continue')
            self.__Engine._VisItBlockingComm = False
            self.__Engine._VisItControlStepping = False
            self.__Engine._UpdateVisItGUI = False

        else:
            pError('Unrecgnozied command %s'%(qcmd))

        return

##############################################################################
class MonitorBehavior(object):
    """ Implements monitor behavior """
    #-------------------------------------------------------------------------
    def __init__(self, engine):
        """ """
        self.__Engine = engine
        return

    #-------------------------------------------------------------------------
    def Initialize(self):
        """
        Initialize the object for monitor mode
        """
        self.__Engine.Disconnect()
        self.__Engine._UpdateVisItGUI = False
        self.__Engine._VisItControlStepping = False
        self.__Engine._VisItBlockingComm = False
        self.__Engine._Synchronous = 0
        self.__Engine._CommandQueue = []
        return

    #-------------------------------------------------------------------------
    def Connect(self):
        """
        Respond to connect in monitor mode
        """
        self.__Engine._UpdateVisItGUI = True
        self.__Engine._VisItControlStepping = True
        self.__Engine._VisItBlockingComm = True
        self.__Engine._Synchronous = 0
        self.__Engine._CommandQueue = []
        return

    #-------------------------------------------------------------------------
    def Update(self):
        """
        Respond to request for simulation step in monitor mode
        """
        # interactive update
        if self.__Engine._Connected:
            # render plot for this sim time step
            if self.__Engine._UpdateVisItGUI:
                # trigger render of plot defined in the GUI
                simV2.VisItUpdatePlots()

        # advance the simulation if auto stepping
        if not self.__Engine._VisItControlStepping:

            # not already rendering
            # check if the sim is finished
            if self.__Engine._Simulation.Continue():
                # sim is not finished
                # take a step
                self.__Engine.StepSimulation()
                self.__Engine.ProbeMemory()
            else:
                # sim is finished
                # send command to end the simulation
                self.__Engine.RequestShutdown()

        return

    #-------------------------------------------------------------------------
    def Process(self, qcmd):
        """
        process monitor control commands.
        """
        if qcmd == 'pause':
            pStatus('WarpVisItEngine pause')
            self.__Engine._VisItBlockingComm = True
            self.__Engine._VisItControlStepping = True
            self.__Engine._UpdateVisItGUI = True
            self.Update()

        elif qcmd == 'step':
            pStatus('WarpVisItEngine step')
            self.__Engine._VisItBlockingComm = True
            self.__Engine._VisItControlStepping = True
            self.__Engine._UpdateVisItGUI = True
            self.__Engine.StepSimulation()
            self.Update()

        elif qcmd == 'run':
            pStatus('WarpVisItEngine run')
            self.__Engine._VisItBlockingComm = False
            self.__Engine._VisItControlStepping = False
            self.__Engine._UpdateVisItGUI = True

        elif qcmd == 'continue':
            pStatus('WarpVisItEngine continue')
            self.__Engine._VisItBlockingComm = False
            self.__Engine._VisItControlStepping = False
            self.__Engine._UpdateVisItGUI = False

        elif qcmd == 'disconnect':
            pStatus('WarpVisItEngine disconnect')
            self.__Engine._VisItBlockingComm = False
            self.__Engine._VisItControlStepping = False
            self.__Engine._UpdateVisItGUI = False
            self.__Engine.Disconnect()

        else:
            pError('Unrecgnozied command %s'%(qcmd))
        return

##############################################################################
class PromptBehavior(object):
    """ Implements command prompt behavior """
    #-------------------------------------------------------------------------
    def __init__(self, engine):
        """ """
        self.__Engine = engine
        self.__Globals = {}
        # for prompt mode
        # TODO -- not working when run in MPI
        import readline
        readline.parse_and_bind("tab: complete")
        return

    #-------------------------------------------------------------------------
    def Initialize(self):
        """
        Initialize the object for monitor mode
        """
        self.__Engine.Disconnect()
        self.__Engine._UpdateVisItGUI = False
        self.__Engine._VisItControlStepping = False
        self.__Engine._VisItBlockingComm = False
        self.__Engine._Synchronous = 0
        self.__Engine._CommandQueue = []
        if self.__Engine._Simulation is not None:
            self.__Globals = self.__Engine._Simulation.GetNamespace()
        return

    #-------------------------------------------------------------------------
    def Connect(self):
        """
        Respond to connect in monitor mode
        """
        self.__Engine._UpdateVisItGUI = True
        self.__Engine._VisItControlStepping = True
        self.__Engine._VisItBlockingComm = True
        self.__Engine._Synchronous = 0
        self.__Engine._CommandQueue = []
        return

    #-------------------------------------------------------------------------
    def Update(self):
        """
        Respond to request for simulation step in monitor mode
        """
        masterRank = lambda r: r==0
        if self.__Engine._Connected:
            # interactive update
            # render plot for this sim time step
            if self.__Engine._UpdateVisItGUI:
                # trigger render of plot defined in the GUI
                simV2.VisItUpdatePlots()

            # advance the simulation if auto stepping
            if not self.__Engine._VisItControlStepping:
                self.__Engine.StepSimulation()

        else:
            # get input from user
            cmd=''
            if masterRank(self.__Engine._CommRank):
                try:
                    cmd = raw_input('> ')
                except:
                    einf = sys.exc_info()
                    traceback.print_exception(einf[0], einf[1], None)
            cmd = parallel.broadcast(cmd, root=0)
            self.__Engine.ProcessCommand(cmd)

        return

    #-------------------------------------------------------------------------
    def Process(self, qcmd):
        """
        process prompt control commands.
        """
        masterRank = lambda r: r==0
        if self.__Engine._Connected:
            # hanlde commands from the GUI
            if qcmd == 'pause':
                pStatus('WarpVisItEngine pause')
                self.__Engine._VisItBlockingComm = True
                self.__Engine._VisItControlStepping = True
                self.__Engine._UpdateVisItGUI = True
                self.Update()

            elif qcmd == 'step':
                pStatus('WarpVisItEngine step')
                self.__Engine._VisItBlockingComm = True
                self.__Engine._VisItControlStepping = True
                self.__Engine._UpdateVisItGUI = True
                self.__Engine.StepSimulation()
                self.Update()

            elif qcmd == 'run':
                pStatus('WarpVisItEngine run')
                self.__Engine._VisItBlockingComm = False
                self.__Engine._VisItControlStepping = False
                self.__Engine._UpdateVisItGUI = True

            elif qcmd == 'continue':
                pStatus('WarpVisItEngine continue')
                self.__Engine._VisItBlockingComm = False
                self.__Engine._VisItControlStepping = False
                self.__Engine._UpdateVisItGUI = False

            elif qcmd == 'end':
                pStatus('WarpVisItEngine disconnect')
                self.__Engine._VisItBlockingComm = False
                self.__Engine._VisItControlStepping = False
                self.__Engine._UpdateVisItGUI = False
                self.__Engine.Disconnect()

        else:
            # handle user prompt
            if (qcmd == 'quit'):
                self.__Engine.RequestShutdown()

            elif qcmd == 'connect':
                pStatus('WarpVisItEngine connect')
                # put the engine in interact mode
                self.__Engine._UpdateVisItGUI = False
                self.__Engine._VisItControlStepping = True
                self.__Engine._VisItBlockingComm = True
                self.__Engine._Synchronous = 0
                self.__Engine._CommandQueue = []
            else:
                try:
                    exec qcmd in self.__Globals
                except:
                    if masterRank(self.__Engine._CommRank):
                        einf = sys.exc_info()
                        traceback.print_exception(einf[0], einf[1], None)
                    pass

        return

##############################################################################
class BatchBehavior(object):
    """ Implements batch mode behvaior """
    #-------------------------------------------------------------------------
    def __init__(self, engine):
        """ """
        self.__Engine = engine
        return

    #-------------------------------------------------------------------------
    def Initialize(self):
        """
        Initialize the object for batch mode
        """
        self.__Engine.Disconnect()
        self.__Engine._UpdateVisItGUI = False
        self.__Engine._VisItControlStepping = False
        self.__Engine._VisItBlockingComm = True
        self.__Engine._Synchronous = 0
        self.__Engine._CommandQueue = []
        return

    #-------------------------------------------------------------------------
    def Connect(self):
        """
        Respond to connect in batch mode
        """
        self.__Engine._UpdateVisItGUI = False
        self.__Engine._VisItControlStepping = False
        self.__Engine._VisItBlockingComm = True
        self.__Engine._Synchronous = 0
        self.__Engine._CommandQueue = []
        self.__Engine.ProbeMemory(True)
        return

    #-------------------------------------------------------------------------
    def Update(self):
        """
        Respond to request for simulation step in batch mode
        """
        if not self.__Engine.GetSynchronous():
            # not already rendering
            # render plots for this sim time step
            if self.__Engine._Connected:
                self.__Engine.PushRenderScripts()

            # not already rendering
            # check if the sim is finished
            if self.__Engine._Simulation.Continue():
                # sim is not finished
                # send command to take a step
                self.__Engine.ProcessCommand('adv')

            else:
                # sim is finished
                # send command to end the simulation
                self.__Engine.ProcessCommand('end')
        return

    #-------------------------------------------------------------------------
    def Process(self, qcmd):
        """
        process batch control commands.

        adv : take a step and intiate the next step
        end : shut everything down
        """
        if qcmd == 'adv':
            self.__Engine.StepSimulation()
            self.__Engine.ProbeMemory()
            self.Update()

        elif (qcmd == 'end'):
            self.__Engine.RequestShutdown()

        else:
            pError('Unrecgnozied command %s'%(qcmd))
        return



#-----------------------------------------------------------------------------
def broadcastInt(val, sender=0) :
    """Integer broadcast callback for libsim"""
    result = parallel.broadcast(val, root=sender)
    if __debug__: pDebug('broadcastInt %i'%(val))
    return result

#-----------------------------------------------------------------------------
def broadcastString(val, n=0, sender=0) :
    """String broadcast callback for libsim"""
    result = parallel.broadcast(val, root=sender)
    if __debug__: pDebug('broadcastString %s'%(result))
    return result

#-----------------------------------------------------------------------------
def broadcastSlaveCommand(command=0):
    """Helper function for ProcessCommands"""
    result = parallel.broadcast(command, root=0)
    if __debug__: pDebug('broadcastSlaveCommand %i'%(result))
    return result

#-----------------------------------------------------------------------------
def slaveProcessCallback():
    """Callback involved in command communication"""
    if __debug__: pDebug('slaveProcessCallback')
    processCommand = int(0)
    broadcastSlaveCommand(processCommand)

#-----------------------------------------------------------------------------
def commandCallback(cmd, args, engine):
    """
    Callback function used to process simulation
    control commands.
    """
    engine.ProcessCommand(cmd)
    return
