from WarpVisItUtil import pError
from WarpVisItUtil import pDebug
import parallel
import multiprocessing
import visit
from visit import visit
import os
import socket

#-----------------------------------------------------------------------------
class WarpVisItCLI:
    """
    Helper class which for non-interactive runs forks a CLI
    which in turn launches a Viewer. WarpVisIt CLI manages
    the forked processes.
    """
    #-------------------------------------------------------------------------
    def __init__(self):
        """ """
        self.__CommRank = parallel.get_rank()
        self.__CLIRank = 0
        self.__CLIProc = None
        self.__ListenSock = None
        self.__CLISock = None
        self.__SimFile = None
        self.__RenderCallback = None
        self.__ShouldRenderCallback = None
        self.__CLIHost = 'localhost'
        self.__CLIPort = 49494

    #-------------------------------------------------------------------------
    def SetSimFile(self, fileName):
        """
        Set the .sim2 file name.
        """
        self.__SimFile = fileName

    #-------------------------------------------------------------------------
    def SetRenderCallback(self, f):
        """
        Set the callback used to render results. If one is provided
        then the view will run non-interactively.
        """
        self.__RenderCallback = f

    #-------------------------------------------------------------------------
    def SetShouldRenderCallback(self, f):
        """
        Set the callback used to determine when to render results.
        """
        self.__ShouldRenderCallback = f

    #-------------------------------------------------------------------------
    def GetSimFile(self):
        """Return the .sim2 file name"""
        return self.__SimFile

    #-------------------------------------------------------------------------
    def EventLoop(self):
        """
        """
        pDebug('WarpVisItCLI::EventLoop')
        if ((self.__RenderCallback is not None)
                and (self.__CLIRank == self.__CommRank)):
            # ask the simulation if we should render
            # this step. transfer the result to the
            # CLI which then in turn forwards it to
            # the viewer.
            if self.__ShouldRenderCallback():
                self.__CLISock.sendall('render')

            else:
                self.__CLISock.sendall('skip')

        return True

    #-------------------------------------------------------------------------
    def Initialize(self, args=['-nowin']):
        """
        Open the VisIt viewer, connect to the simulation and then execute
        the visualization defined by the visFunction. Any given arguments
        are passed on to the viewer in the form of command line arguments.
        """
        pDebug('WarpVisItCLI::Initialize')
        if ((self.__RenderCallback is not None)
                and (self.__CLIRank == self.__CommRank)):
            # open a communication channel for the newly spawned viewer
            self._ListenForCLI()

            # launch the CLI which in turn launches the viewer.
            # it's in a separate process, this is VisIt's design.
            self.__CLIProc = multiprocessing.Process(
                  target=CLIMain,
                  args=(self.__SimFile,
                  self.__RenderCallback,
                  self.__CLIHost,
                  self.__CLIPort,
                  args))

            # daemon means automatic cleanup of child processes
            self.__CLIProc.daemon = True
            self.__CLIProc.start()

            self.__CLISock, addr = self.__ListenSock.accept()
            pDebug('Connected on %s'%(str(addr)))

        return True

    def Finalize(self):
        """
        Shutdown clean up etc...
        """
        pDebug('WarpVisItCLI::Finalize')
        if self.__CLIRank == self.__CommRank:

            self.__CLISock.sendall('quit')
            self.__CLISock.close()
            self.__CLISock = None

            self.__ListenSock.close()
            self.__ListenSock = None

            #self.__CLIProc.terminate()
            #self.__CLIProc = None

        return

    #-------------------------------------------------------------------------
    def _ListenForCLI(self):
        """
        Create a socket connection on host:port
        """
        pDebug('WarpVisItCLI::Connect')
        for candidate in socket.getaddrinfo(
                self.__CLIHost,
                self.__CLIPort,
                socket.AF_UNSPEC,
                socket.SOCK_STREAM,
                0,
                socket.AI_PASSIVE):

            pDebug('trying %s'%(str(candidate)))

            af, socktype, proto, canonname, sa = candidate

            try:
                self.__ListenSock = socket.socket(af, socktype, proto)

            except socket.error as msg:
                pDebug('failed to open %s'%(str(candidate)))
                self.__ListenSock = None
                continue

            try:
                self.__ListenSock.bind(sa)
                self.__ListenSock.listen(1)

            except socket.error as msg:
                pDebug('failed to bind %s'(str(candidate)))
                self.__ListenSock.close()
                self.__ListenSock = None
                continue

            break

        if self.__ListenSock is None:
            pError('Failed to create viewer socket on %s:%d'%(host,port))
            return False

        pDebug('sucess')
        return True

#-----------------------------------------------------------------------------
def CLIMain(simFile, renderCallback, host, port, args=[]):
    """
    This is the viewer main loop. The viewer resides in a separate process.
    """
    from visit import visit
    import socket

    # CLI establish connection with the engine
    engineSock = None
    for candidate in socket.getaddrinfo(
            host,
            port,
            socket.AF_UNSPEC,
            socket.SOCK_STREAM):

        print 'CLI main trying %s'%(str(candidate))

        af, socktype, proto, canonname, sa = candidate

        try:
            engineSock = socket.socket(af, socktype, proto)
        except socket.error as msg:
            engineSock = None
            continue

        try:
            engineSock.connect(sa)
        except socket.error as msg:
            engineSock.close()
            engineSock = None
            continue

        break

    if engineSock is None:
        raise RuntimeError(
        'CLI failed to connect to the engine at %s:%d'%(host,port))

    print 'CLI connected to the engine at %s:%d'%(host,port)

    # Were actually in the CLI, here is where we start the viewer
    # process
    for arg in args:
        visit.AddArgument(arg)
    visit.Launch()
    #visit.ShowAllWindows()
    visit.OpenDatabase(simFile)

    # We in CLI drive both viewer and engine
    # but we want to coordinate with the engine
    # because the engine has access to the simulation
    # state. The simulation can then tell us when to
    # render something or even when to terminate.
    render = lambda cmd: cmd=='render'
    skip = lambda cmd: cmd=='skip'
    quit = lambda cmd: cmd=='quit'
    while 1:
        # in non-interactive runs the game is to drive
        # the simulation. When somehtin interesting occurs
        # the simulation sends us a render command.
        visit.SendSimulationCommand('localhost', simFile, 'step')
        cmd = engineSock.recv(1024)
        print 'CLI command %s recvd'%(cmd)
        if skip(cmd):
            continue

        elif render(cmd):
            renderCallback()
            continue

        elif quit(cmd):
            break

    # CLI is finished
    print 'CLI finished'
    exit(0)
