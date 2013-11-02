import os
import simV2
import WarpVisItSimV2Db
from WarpVisItUtil import pError
from WarpVisItUtil import pDebug
from WarpVisItUtil import VisItEnv
from WarpVisItCLI import WarpVisItCLI
import parallel
import warp

#-----------------------------------------------------------------------------
def broadcastInt(val, sender=0) :
    """Integer broadcast callback for libsim"""
    result = parallel.broadcast(val, root=sender)
    pDebug('broadcastInt %i'%(val))
    return result

#-----------------------------------------------------------------------------
def broadcastString(val, n=0, sender=0) :
    """String broadcast callback for libsim"""
    result = parallel.broadcast(val, root=sender)
    pDebug('broadcastString %s'%(result))
    return result

#-----------------------------------------------------------------------------
def broadcastSlaveCommand(command=0):
    """Helper function for ProcessCommands"""
    result = parallel.broadcast(command, root=0)
    pDebug('broadcastSlaveCommand %i'%(result))
    return result

#-----------------------------------------------------------------------------
def slaveProcessCallback():
    """Callback involved in command communication"""
    pDebug('slaveProcessCallback')
    processCommand = int(0)
    broadcastSlaveCommand(processCommand)

#-----------------------------------------------------------------------------
def commandCallback(cmd, args, userData):
    """
    Callback function used to process user-defined simulation
    control commands issues by the VisIt viewer.
    """
    pDebug('commandCallback %s'%(cmd))
    if cmd == "halt":
        userData.SetVisItControl(True)
        userData.SetVisItUpdates(True)
        simV2.VisItUpdatePlots()
        return True

    elif cmd == "step":
        userData.SetVisItControl(True)
        userData.SetVisItUpdates(True)
        cont = userData.Update()
        return cont

    elif cmd == "run":
        userData.SetVisItControl(False)
        userData.SetVisItUpdates(True)
        return True

    elif cmd == "run_without_update":
        userData.SetVisItControl(False)
        userData.SetVisItUpdates(False)
        return True

    elif cmd == "end":
        userData.Abort()
        return True

    pError('Unrecgnozied command %s'%(cmd))
    return False







#-----------------------------------------------------------------------------
class WarpVisItEngine:
    """
    The WarpVisItEngine manages engine state and exposes
    engine API specifically to the Warp simulation.
    """
    COMMAND_PROCESS = 0
    COMMAND_SUCCESS = 1
    COMMAND_FAILURE = 2

    #-------------------------------------------------------------------------
    def __init__(self):
        """ """
        self.__CommRank = parallel.get_rank()
        self.__CommSize = parallel.number_of_PE()
        self.__Env = VisItEnv()
        self.__SimFile = ''
        self.__HasTrace = False
        self.__MeshCallback = None
        self.__MetaDataCallback = None
        self.__VariableCallback = None
        self.__DomainListCallback = None
        self.__TimeCallback = None
        self.__IterationCallback = None
        self.__StepCallback = None
        self.__StepInterval = 1
        self.__ContinueCallback = None
        self.__ActiveRenderScriptsCallback = None
        self.__RenderScripts = {}
        self.__VisItBlockingComm = True
        self.__VisItUpdates = False
        self.__VisItControl = True
        self.__VisItAbort = False
        self.__Interactive = True

        simV2.VisItSetDirectory(self.__Env.GetRoot())
        simV2.VisItSetupEnvironment()
        simV2.VisItSetParallelRank(self.__CommRank)
        simV2.VisItSetParallel(self.__CommSize > 1)
        simV2.VisItSetBroadcastIntFunction(broadcastInt)
        simV2.VisItSetBroadcastStringFunction(broadcastString)

        self.SetMetaDataCallback(WarpVisItSimV2Db.getMetaData)
        self.SetMeshCallback(WarpVisItSimV2Db.getMesh)
        self.SetVariableCallback(WarpVisItSimV2Db.getVar)
        self.SetDomainListCallback(WarpVisItSimV2Db.getDomains)
        self.SetTimeCallback(WarpVisItSimV2Db.getTime)
        self.SetIterationCallback(WarpVisItSimV2Db.getIteration)

    #-------------------------------------------------------------------------
    def SetMeshCallback(self, f):
        """Set the callback that creates meshes"""
        self.__MeshCallback = f

    #-------------------------------------------------------------------------
    def SetMetaDataCallback(self, f):
        """Set the callback that described datasets"""
        self.__MetaDataCallback = f

    #-------------------------------------------------------------------------
    def SetVariableCallback(self, f):
        """Set the callback that creates variables"""
        self.__VariableCallback = f

    #-------------------------------------------------------------------------
    def SetDomainListCallback(self, f):
        """Set the callback that describes domains"""
        self.__DomainListCallback = f

    #-------------------------------------------------------------------------
    def SetContinueCallback(self, f):
        """
        Set the callback that returns true while the simulation should run
        """
        self.__ContinueCallback = f

    #-------------------------------------------------------------------------
    def SetStepCallback(self, f):
        """Set the callback used to advance the simulaition"""
        self.__StepCallback = f

    #-------------------------------------------------------------------------
    def SetActiveRenderScriptsCallback(self, f):
        """
        Set the callback used to query the simulation to check if
        rendering is needed at this time step
        """
        self.__ActiveRenderScriptsCallback = f

    #-------------------------------------------------------------------------
    def SetRenderScripts(self, scripts):
        """
        Set the dictionary of rendering scripts.
        """
        self.__RenderScripts = scripts

    #-------------------------------------------------------------------------
    def SetTimeCallback(self, f):
        """Set the callback used to get the simulation time"""
        self.__SimulationTimeCallback = f

    #-------------------------------------------------------------------------
    def SetIterationCallback(self, f):
        """Set the callback used to get the simulation iteration"""
        self.__SimulationIterationCallback = f

    #-------------------------------------------------------------------------
    def SetStepInterval(self, n):
        """Set the number of simuylation steps to take per VisIt step"""
        self.__StepInterval = int(n)

    #-------------------------------------------------------------------------
    def GetStepInterval(self):
        """Set the number of simuylation steps to take per VisIt step"""
        return self.__StepInterval

    #-------------------------------------------------------------------------
    def SetVisItBlockingComm(self, v):
        """Enable blocking on VisIt commands in the event loop"""
        self.__VisItBlockingComm = bool(v)

    #-------------------------------------------------------------------------
    def SetVisItUpdates(self, v):
        """Enable VisIt rendering"""
        self.__VisItUpdates = bool(v)

    #-------------------------------------------------------------------------
    def SetVisItControl(self, v):
        """Select between interactive/non-interactive simulation stepping"""
        self.__VisItControl = bool(v)
        self.__VisItBlockingComm = bool(v)

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
        self.__VisItUpdates = True if self.__Interactive else False

    #-------------------------------------------------------------------------
    def Initalize(self, simName='Warp', simComment=None, simPath=None,
        traceFile=None, visitOpts=''):
        """
        Perform the initial setup of VisIt to include the libsim module
        with the simulation.
        """
        pDebug('WarpVisItEngine::Initialize')
        if not simPath:
            simPath = os.getcwd()

        if not self.__SimFile:
            self.__SimFile = os.path.join(simPath, simName + '.sim2')

        pDebug('sim2 file = %s'%(self.__SimFile))

        simV2.VisItSetDirectory(self.__Env.GetRoot())

        if traceFile:
            VisItOpenTraceFile(visitTraceFile)
            self.__HasTrace = True

        if visitOpts:
            VisItSetOptions(visitOpts);

        if self.__CommRank == 0:
            if not simV2.VisItInitializeSocketAndDumpSimFile(
                  simName,
                  simComment,
                  simPath,
                  None,
                  None,
                  self.__SimFile):
                pError('VisItInitializeSocketAndDumpSimFile failed')

        if not self.__Interactive:
            pDebug('intializing CLI')
            cli = WarpVisItCLI()
            cli.SetSimFile(self.GetSimFile())
            cli.Initialize()

        return

    #-------------------------------------------------------------------------
    def Finalize(self):
        """Finalize VisIt libsim"""
        pDebug('WarpVisItEngine::Finalize')
        if self.__HasTrace:
            simV2.VisItCloseTraceFile()

    #-------------------------------------------------------------------------
    def ProcessCommands(self):
        """Process commands from the VisIt viewer on all MPI processes"""
        masterRank = lambda r: r==0

        pDebug('WarpVisItEngine::ProcessCommands')

        if masterRank(self.__CommRank):
            if simV2.VisItProcessEngineCommand():
                pDebug('master success')
                broadcastSlaveCommand(self.COMMAND_SUCCESS)
                return True

            else:
                pDebug('master discconnect')
                broadcastSlaveCommand(self.COMMAND_FAILURE)
                return False

        else:
            while True:
                command = broadcastSlaveCommand()
                pDebug('slave command %i'%(command))

                if command == self.COMMAND_PROCESS:
                    pDebug('slave process command')
                    simV2.VisItProcessEngineCommand()

                elif command == self.COMMAND_SUCCESS:
                    pDebug('slave success')
                    return True

                elif command == self.COMMAND_FAILURE:
                    pDebug('slave disconnect')
                    return False

                else:
                    pError('Bad command %i'%(command))

        return False

    #-------------------------------------------------------------------------
    def EventLoop(self):
        """Run the VisIt controlled event loop. Return False on error"""
        doError = lambda e: e<0
        doUpdate = lambda e: e==0
        doConnect = lambda e: e==1
        doCommand = lambda e: e==2
        masterRank = lambda r: r==0

        pDebug('WarpVisItEngine::EventLoop')

        while 1:
            if self.__VisItAbort:
                pDebug('abort')
                return False

            if masterRank(self.__CommRank):
                event = simV2.VisItDetectInput(self.__VisItBlockingComm, -1)
                parallel.broadcast(event)
            else:
                event = parallel.broadcast(0)

            pDebug('event=%d'%(event))

            if doError(event):
                pError('VisItDetectInput detected an error')
                return False

            elif doUpdate(event):
                pDebug('update')
                if self.Update():
                    continue
                return True

            elif doConnect(event):
                pDebug('connect')
                if simV2.VisItAttemptToCompleteConnection() == 1:
                    simV2.VisItSetGetMetaData(self.__MetaDataCallback, self)
                    simV2.VisItSetGetMesh(self.__MeshCallback, self)
                    simV2.VisItSetGetVariable(self.__VariableCallback, self)
                    simV2.VisItSetGetDomainList(self.__DomainListCallback, self)
                    simV2.VisItSetCommandCallback(commandCallback, self)
                    simV2.VisItSetSlaveProcessCallback(slaveProcessCallback)
                    # send ourself a command to kick off the first update
                    # after we update we'll continue to send our self a
                    # a command to kick of the next update.
                    if not self.__Interactive:
                        pDebug('initializing stepping')
                        simV2.VisItExecuteCommand(
                            ("from visit import visit\n"
                             "visit.SendSimulationCommand('localhost', '%s', 'step')\n")%(
                             self.__SimFile))
                    continue
                else:
                    pError('Connection failed')
                    return False

            elif doCommand(event):
                pDebug('command')
                self.__VisItBlockingComm = True
                if not self.ProcessCommands():
                    simV2.VisItDisconnect()
                    self.__VisItBlockingComm = False
                    pDebug('VisIt viewer disconnected')
                    return False
                continue

            else:
                pDebug('unknown command')

            continue

        return False

    #-------------------------------------------------------------------------
    def Update(self):
        """Advance the simualtion update plots etc..."""
        pDebug('WarpVisItEngine::AdvanceAndUpdate')
        self.StepSimulation()
        self.Render()
        if self.__ContinueCallback():
            if not self.__Interactive:
                # send ourself a command to kick off the next update
                # this keeps everything running without any viewer/CLI
                # code.
                simV2.VisItExecuteCommand(
                    ("from visit import visit\n"
                     "visit.SendSimulationCommand('localhost', '%s', 'step')\n")%(
                    self.__SimFile))
            return True
        return False

    #-------------------------------------------------------------------------
    def Render(self):
        """Render simulation data"""
        pDebug('WarpVisItEngine::Render')

        if self.__Interactive:
            if self.__VisItUpdates:
                simV2.VisItUpdatePlots()

        else:
            for script in self.__ActiveRenderScriptsCallback():
                pDebug('Rendering %s'%(script))
                simV2.VisItExecuteCommand(self.__RenderScripts[script])

    #-------------------------------------------------------------------------
    def StepSimulation(self):
        """Drive the simulation through one or more steps"""
        pDebug('WarpVisItEngine::StepSimulation %d'%(self.__StepInterval))
        i = 0
        while i < self.__StepInterval:
            self.__StepCallback()
            simV2.VisItTimeStepChanged()
            i += 1
