import os
import sys
import argparse
import simV2
import WarpVisItSimV2Db
from WarpVisItUtil import pError, pDebug, pStatus, VisItEnv, getEnvVar
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
        callModeMethod = lambda f : getattr(self, f+self.__Mode)()

        self.__CommRank = parallel.get_rank()
        self.__CommSize = parallel.number_of_PE()
        self.__Env = VisItEnv()
        self.__SimFile = getEnvVar('WARPVISIT_SIM_FILE',str,'WarpVisIt.sim2')
        self.__TraceFile = getEnvVar('WARPVISIT_TRACE_FILE',bool,False)
        self.__Simulation = None
        self.__ShutdownRequested = False
        self.__Connected = False
        self.__Synchronous = 0
        self.__Mode = 'Batch'
        self.__UpdateVisItGUI = True
        self.__VisItControlStepping = True
        self.__VisItBlockingComm = True
        self.__CommandQueue = []
        self.__EngineOpts = None
        self.__StepCount = 0
        self.__ProbeMem = 0

        # parse command line args
        interact = getEnvVar('WARPVISIT_MODE_INTERACT',bool,False)
        monitor = getEnvVar('WARPVISIT_MODE_MONITOR',bool,False)
        batch = getEnvVar('WARPVISIT_MODE_BATCH',bool,False)

        ap = argparse.ArgumentParser(usage=argparse.SUPPRESS,prog='WarpVisItEngine',add_help=False)
        ap.add_argument('--sim-file',type=str,default=self.__SimFile)
        ap.add_argument('--trace-file',default=self.__TraceFile, action='store_true')
        ap.add_argument('--interact',default=interact, action='store_true')
        ap.add_argument('--monitor',default=monitor, action='store_true')
        ap.add_argument('--batch',default=batch, action='store_true')
        ap.add_argument('--engine-opts',type=str,default=None)
        ap.add_argument('--probe-mem',type=int,default=self.__ProbeMem)
        opts = vars(ap.parse_known_args(args)[0])
        self.__SimFile = os.path.abspath(opts['sim_file'])
        self.__TraceFile = opts['trace_file']
        interact = opts['interact']
        monitor = opts['monitor']
        batch = opts['batch']
        self.__EngineOpts = opts['engine_opts']
        self.__ProbeMem = opts['probe_mem']

        if sum((interact,monitor,batch))>1:
            raise RuntimeError('--interact, --monitor, --batch are mutually exclusive')
        if interact:
            self.__Mode = 'Interact'
        elif monitor:
            self.__Mode = 'Monitor'
        else:
            self.__Mode = 'Batch'

        callModeMethod('Initialize')

        # initialize libsim
        simV2.VisItSetDirectory(self.__Env.GetRoot())
        simV2.VisItSetupEnvironment()
        simV2.VisItSetParallelRank(self.__CommRank)
        simV2.VisItSetParallel(self.__CommSize > 1)
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
        self.__Simulation = simulation

    #-------------------------------------------------------------------------
    def RequestShutdown(self):
        """Request shutdown"""
        self.__VisItControlStepping = False
        self.__VisItBlockingComm = False
        self.__UpdateVisItGUI = False
        self.__ShutdownRequested = True

    #-------------------------------------------------------------------------
    def SetVisItControlStepping(self, v):
        """Select between interactive/non-interactive simulation stepping"""
        self.__VisItControlStepping = bool(v)
        self.__VisItBlockingComm = bool(v)

    #-------------------------------------------------------------------------
    def GetVisItControlStepping(self, v):
        """Select between interactive/non-interactive simulation stepping"""
        return self.__VisItControlStepping

    #-------------------------------------------------------------------------
    def SetSynchronous(self, v):
        """
        Tracks wether we are currently rendering. we need to prevent simulation
        advances while rendering, else we may mix data from multiple time steps
        in the same plots.
        """
        if bool(v):
            self.__Synchronous += 1
        else:
            self.__Synchronous -= 1
        if self.__Synchronous < 0:
            self.__Synchronous = 0

    #-------------------------------------------------------------------------
    def GetSynchronous(self):
        """True when VisIt is actively rendering"""
        return self.__Synchronous > 0

    #-------------------------------------------------------------------------
    def SetSimFile(self, fileName):
        """
        Set the .sim2 file name. If one is not set the default is
        '/cwd/Warp.sim2'
        """
        self.__SimFile = fileName

    #-------------------------------------------------------------------------
    def GetSimFile(self):
        """Return the .sim2 file name"""
        return self.__SimFile

    #-------------------------------------------------------------------------
    def SetInteractive(self, interactive):
        """Set true if the simulation should connect to VisIt"""
        self.__Interactive = bool(interactive)
        self.__UpdateVisItGUI = bool(interactive)

    #-------------------------------------------------------------------------
    def Initalize(self, engineOpts=''):
        """
        Perform the initial setup of VisIt to include the libsim module
        with the simulation.
        """
        if __debug__: pDebug('WarpVisItEngine::Initialize')

        simV2.VisItSetDirectory(self.__Env.GetRoot())

        if self.__TraceFile:
            simV2.VisItOpenTraceFile('%d.trace'%(self.__CommRank))

        if self.__EngineOpts:
            if __debug__: pDebug('Using options %s'%(self.__EngineOpts))
            simV2.VisItSetOptions(self.__EngineOpts);

        if self.__CommRank == 0:
            if not simV2.VisItInitializeSocketAndDumpSimFile(
                  'WarpVisIt', None, os.getcwd(), None, None, self.__SimFile):
                pError('VisIt initialization failed')
                return False

        return True

    #-------------------------------------------------------------------------
    def Finalize(self):
        """
        Finalize VisIt libsim, cleanup, etc
        """
        if self.__TraceFile:
            simV2.VisItCloseTraceFile()
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
        onErrorContinue = lambda : self.__Mode=='Monitor'
        callModeMethod = lambda f : getattr(self, f+self.__Mode)()

        if __debug__: pDebug('WarpVisItEngine::EventLoop')

        while 1:
            # test or wait for incomming communication.
            if masterRank(self.__CommRank):
                event = simV2.VisItDetectInput(self.__VisItBlockingComm, -1)
                parallel.broadcast(event)
            else:
                event = parallel.broadcast(0)

            # process incomming communication
            if __debug__: pDebug('event=%d'%(event))

            # internal comm error
            if doError(event):
                # the VisIt viewer disconnected.
                pError('An communication error was detected')
                callModeMethod('Initialize')
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
                callModeMethod('Update')

            # internal connection
            elif doConnect(event):
                if simV2.VisItAttemptToCompleteConnection() == 1:
                    pStatus('WarpVisItEngine connected')
                    self.ConnectLibsim()
                    callModeMethod('Connect')
                    callModeMethod('Update')
                else:
                    pError('Connection failed')
                    callModeMethod('Initialize')
                    if onErrorContinue():
                        continue
                    else:
                        return False

            # internal comm
            elif doInternal(event):
                if __debug__: pDebug('internal')
                if not self.CommunicateLibsim():
                    pError('Error while processing internal commands')
                    callModeMethod('Initialize')
                    if onErrorContinue():
                        continue
                    else:
                        return False

            # bad value
            else:
                pError('unknown command %d'%(event))

            # check for asynchronous shutdown request
            if self.__ShutdownRequested:
                pStatus('WarpVisItEngine exiting the event loop')
                simV2.VisItDisconnect()
                self.__Connected = False
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
        self.__Connected = True
        self.__VisItBlockingComm = True
        return

    #-------------------------------------------------------------------------
    def CommunicateLibsim(self):
        """ProcessCommand internal VisIt commands."""
        masterRank = lambda r: r==0
        COMMAND_PROCESS = 0
        COMMAND_SUCCESS = 1
        COMMAND_FAILURE = 2

        if __debug__: pDebug('WarpVisItEngine::ProcessCommands')

        if masterRank(self.__CommRank):
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

        adv       : take a step and intiate the next step
        step      : take one simulation step and wait for more commands
        run       : make sequential steps and poll for commands
        continue  : like run but do not update plots
        pause     : finish the current step and wait for more commands
        end       : shut everything down
        endSyn    : internal use only, for synchronziation between sim and vis

        """
        callModeMethod = lambda f : getattr(self, f+self.__Mode)()

        if __debug__: pDebug('WarpVisItEngine::ProcessCommand %s'%(cmd))

        self.__CommandQueue.append(cmd)

        # queue commands while in synchronous mode
        if self.GetSynchronous():
            if cmd == 'endSyn':
                self.SetSynchronous(False)
                self.__CommandQueue.pop()
            else:
                if __debug__: pDebug('defered %s'%(cmd))
                return

        # now in asynchronous mode
        # process queued commands
        ncmds = len(self.__CommandQueue)
        i = 0
        while (i < ncmds):
            qcmd = self.__CommandQueue.pop()
            i += 1

            if __debug__: pDebug('processes %s'%(qcmd))

            if qcmd == 'adv':
                self.StepSimulation()
                self.ProbeMemory()
                callModeMethod('Update')

            elif (qcmd == 'end'):
                self.RequestShutdown()

            elif qcmd == 'pause':
                pStatus('WarpVisItEngine pause')
                self.__VisItBlockingComm = True
                self.__VisItControlStepping = True
                self.__UpdateVisItGUI = True
                callModeMethod('Update')

            elif qcmd == 'step':
                pStatus('WarpVisItEngine step')
                self.__VisItBlockingComm = True
                self.__VisItControlStepping = True
                self.__UpdateVisItGUI = True
                self.StepSimulation()
                callModeMethod('Update')

            elif qcmd == 'run':
                pStatus('WarpVisItEngine run')
                self.__VisItBlockingComm = False
                self.__VisItControlStepping = False
                self.__UpdateVisItGUI = True

            elif qcmd == 'continue':
                pStatus('WarpVisItEngine continue')
                self.__VisItBlockingComm = False
                self.__VisItControlStepping = False
                self.__UpdateVisItGUI = False

            elif qcmd == 'disconnect':
                pStatus('WarpVisItEngine disconnect')
                self.__VisItBlockingComm = False
                self.__VisItControlStepping = False
                self.__UpdateVisItGUI = False
                self.Disconnect()

            else:
                pError('Unrecgnozied command %s'%(qcmd))

        return

    #-------------------------------------------------------------------------
    def PushRenderScripts(self):
        """Render simulation data"""
        if __debug__: pDebug('WarpVisItEngine::Render')

        activeScripts = self.__Simulation.GetActiveRenderScripts()
        if len(activeScripts):
            scripts = self.__Simulation.GetRenderScripts()
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
                          self.__SimFile)
                simV2.VisItExecuteCommand(source)

        return True

    #-------------------------------------------------------------------------
    def StepSimulation(self):
        """Drive the simulation through one or more steps"""
        if __debug__: pDebug('WarpVisItEngine::StepSimulation')

        self.__Simulation.Advance()
        simV2.VisItTimeStepChanged()
        self.__StepCount += 1
        return True

    #-------------------------------------------------------------------------
    def ProbeMemory(self,force=False):
        """ """
        # report memory consumption
        if force or (self.__ProbeMem > 0) and ((self.__StepCount%self.__ProbeMem) == 0):
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
        if self.__Connected:
            simV2.VisItDisconnect()
            self.__Connected = False
        return


    #-------------------------------------------------------------------------
    def InitializeBatch(self):
        """
        Initialize the object for batch mode
        """
        self.Disconnect()
        self.__UpdateVisItGUI = False
        self.__VisItControlStepping = False
        self.__VisItBlockingComm = True
        self.__Synchronous = 0
        self.__CommandQueue = []
        return

    #-------------------------------------------------------------------------
    def ConnectBatch(self):
        """
        Respond to connect in batch mode
        """
        self.__UpdateVisItGUI = False
        self.__VisItControlStepping = False
        self.__VisItBlockingComm = True
        self.__Synchronous = 0
        self.__CommandQueue = []
        self.ProbeMemory(True)
        return

    #-------------------------------------------------------------------------
    def UpdateBatch(self):
        """
        Respond to request for simulation step in batch mode
        """
        if not self.GetSynchronous():
            # not already rendering
            # render plots for this sim time step
            if self.__Connected:
                self.PushRenderScripts()

            # not already rendering
            # check if the sim is finished
            if self.__Simulation.Continue():
                # sim is not finished
                # send command to take a step
                self.ProcessCommand('adv')

            else:
                # sim is finished
                # send command to end the simulation
                self.ProcessCommand('end')


    #-------------------------------------------------------------------------
    def InitializeMonitor(self):
        """
        Initialize the object for monitor mode
        """
        self.Disconnect()
        self.__UpdateVisItGUI = False
        self.__VisItControlStepping = False
        self.__VisItBlockingComm = False
        self.__Synchronous = 0
        self.__CommandQueue = []
        return

    #-------------------------------------------------------------------------
    def ConnectMonitor(self):
        """
        Respond to connect in monitor mode
        """
        self.__UpdateVisItGUI = True
        self.__VisItControlStepping = True
        self.__VisItBlockingComm = True
        self.__Synchronous = 0
        self.__CommandQueue = []
        return

    #-------------------------------------------------------------------------
    def UpdateMonitor(self):
        """
        Respond to request for simulation step in monitor mode
        """
        # interactive update
        if self.__Connected:
            # render plot for this sim time step
            if self.__UpdateVisItGUI:
                # trigger render of plot defined in the GUI
                simV2.VisItUpdatePlots()

        # advance the simulation if auto stepping
        if not self.__VisItControlStepping:

            # not already rendering
            # check if the sim is finished
            if self.__Simulation.Continue():
                # sim is not finished
                # take a step
                self.StepSimulation()
                self.ProbeMemory()
            else:
                # sim is finished
                # send command to end the simulation
                self.RequestShutdown()

        return

    #-------------------------------------------------------------------------
    def InitializeInteract(self):
        """
        Initialize the object for interact mode
        """
        self.Disconnect()
        self.__UpdateVisItGUI = False
        self.__VisItControlStepping = True
        self.__VisItBlockingComm = True
        self.__Synchronous = 0
        self.__CommandQueue = []
        return

    #-------------------------------------------------------------------------
    def ConnectInteract(self):
        """
        Connect the object for interact mode
        """
        self.__UpdateVisItGUI = True
        self.__VisItControlStepping = True
        self.__VisItBlockingComm = True
        self.__Synchronous = 0
        self.__CommandQueue = []
        return

    #-------------------------------------------------------------------------
    def UpdateInteract(self):
        """
        Respond to request for render in interact mode
        """
        # interactive update
        if self.__Connected:
            # render plot for this sim time step
            if self.__UpdateVisItGUI:
                # trigger render of plot defined in the GUI
                simV2.VisItUpdatePlots()

            # advance the simulation if auto stepping
            if not self.__VisItControlStepping:
                self.StepSimulation()
        return




# TODO these should move into the engine class
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
