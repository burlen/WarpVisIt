import os
import sys
import visit_info
visitlibdir = visit_info.visitdir_lib
sys.path.append(visitlibdir)
from simV2 import *

#Used by warp to wrap mpi function calls
import parallel

#Create a trace file
import tempfile

#Defines
visit_command = { 'VISIT_COMMAND_PROCESS':1 , 'VISIT_COMMAND_SUCCESS':1 , 'VISIT_COMMAND_FAILURE':2 }

def initalize( sim_name="warp", sim_comment="No Simulation Comment", sim_path = os.getcwd() , visitTraceFile = None ) :
    """Perform the initial setup of VisIt to include the libsim module with the simulation
       sim_name : Name for the current simulation. This will appear in the name of the VisIt simulation file.
       sim_comment : A comment string for the current simulation
       sim_path : Path to where the simualtion was started. By default this is assumed to 
                 be the current working directory
       return : string indicating the path to the sim2 file created by libsim    
    """
    #Return parameter
    simFileName = None
    #Set the VisIt home directory
    VisItSetDirectory( visit_info.visitdir )
    options = "-debug 5"
    VisItSetOptions(options);
    #Create the trace file
    deleteTraceFile = False
    if visitTraceFile is  None :
        deleteTraceFile = True
        visitTraceFile = tempfile.mkstemp()[1]
    VisItOpenTraceFile(visitTraceFile)
    #Setyp the visit environment
    VisItSetupEnvironment()
    #Initalize broadcast function for parallel execution
    if parallel.lparallel :
        VisItSetBroadcastIntFunction( broadcast_int )
        VisItSetBroadcastStringFunction( broadcast_string )
        VisItSetParallel( parallel.lparallel)
        VisItSetParallelRank( parallel.get_rank() )
    #Initalize the visit libsim engine
    if parallel.get_rank() == 0 :
        resp = VisItInitializeSocketAndDumpSimFile( sim_name , sim_comment , sim_path , None, None , None )
        #If VisIt libsim was started correctly
        if resp == 1 :
            simFileName = _get_simfilename_from_trace( visitTraceFile )
    #Close and remove the trace file
    if deleteTraceFile : 
        VisItCloseTraceFile()
        os.remove( visitTraceFile )
    return simFileName

def finalize() :
    """Finalize VisIt libsim"""
    VisItCloseTraceFile()

def waitForViewer(visit_datainterface, callback_userdata, timeout=-1) :
    """Blocking function which waits for a VisIt viewer to connect 
    
       Keyword arguments:
       visit_datainterface : python module with the callback functions to be used for VisIt libsim
       callback_userdata  : user data to be communicated to the callback functions
       timeout : interger indicating the mili-seconds to timeout or -1 if the function should not timeout. 
    """
    print "Waiting for visit viewer to connect"
    if timeout <= 0 :
        status = VisItDetectInput( True ,-1 )
    else : 
        status = VisItDetectInputWithTimeout( True , timeout , -1)
    print "VisIt message found by libsim"
    print status
    if status == 1  :
        status = VisItAttemptToCompleteConnection()
        if status==1 : 
            print "VisIt connected"
            connect_callbacks(visit_datainterface, callback_userdata)
            print "VisIt Engine: connected callbacks"
        else :
            print "Connection with VisIt failed"
	return True
    else :
        print "VisIt not connected"
	return False

def connect_callbacks(datainterfaces , callback_userdata) :
    """Connect the necessary VisIt callback functions
    
       The function checks whether a set of callback functions are available in 
       the datainterfaces object provided as input. The following callback functions
       are checked:
       
       get_meta_data : Function used to populate metadata information to VisIt.
                       This function should be implemented in all cases 
       get_mesh : Function used to access mesh data (e.g., particles)
       get_variable : Function used to access variable data, e.g., vx for particles
       control_command_callback : Callback function to handle user-defined control commands issued by the viewer.
    """
    
    callbackList = dir( datainterfaces )
    if 'get_meta_data' in callbackList :
        VisItSetGetMetaData(datainterfaces.get_meta_data, callback_userdata)
        print "Registered GetMetaData callback"
    if 'get_mesh' in callbackList :
        VisItSetGetMesh(datainterfaces.get_mesh, callback_userdata)
        print "Registered GetMesh callback"
    if 'get_variable' in callbackList :
        VisItSetGetVariable(datainterfaces.get_variable, callback_userdata)
        print "Registered GetVariable callback"
    if 'control_command_callback' in callbackList :
        VisItSetCommandCallback(datainterfaces.control_command_callback , callback_userdata )
        print "Registered CommandCallback callback"
    
    #Register callbacks only needed for the parallel case
    if parallel.lparallel : 
        VisItSetSlaveProcessCallback( visit_engine.slave_process_callback )
        if 'get_domain_list' in callbackList :
            VisItSetGetDomainList(datainterfaces.get_domain_list, 0)
            print "Registered GetDomainList callback"


def broadcast_int(ival, sender=0) :
    """Integer broadcast callback for libsim"""
    ret = None
    if parallel.lparallel :
        if parallel.get_rank() == sender :
            ret = mpicom.broadcast(ival)
        else :
            ret = mpicom.broadcast()
    return ret

def broadcast_string(sval, slen, sender=0) :
    """String broadcast callback for libsim"""
    ret = None
    if parallel.lparallel :
        if parallel.get_rank() == sender :
            ret = mpicom.broadcast(ival)
        else : 
            ret = mpicom.broadcast()
    return ret

def broadcast_slave_command( command ) :
    """Helper function for process_visit_command"""
    broadcast_int( command , 0 )

def slave_process_callback() :
    """Callback involved in command communication"""
    command = visit_command['VISIT_COMMAND_PROCESS']
    broadcast_slave_command( command )

def process_visit_command() :
    """Process commands from the VisIt viewer on all MPI processes"""
    if parallel.get_rank() == 0 :
        if VisItProcessEngineCommand() :
            command = visit_command['VISIT_COMMAND_SUCCESS']
            broadcast_slave_command(command)
            return True
        else :
            command = visit_command['VISIT_COMMAND_FAILURE']
            broadcast_slave_command( command )
            return False
    else :
        while True :
            command = 0
            broadcast_slave_command( command ) 
            if command == visit_command['VISIT_COMMAND_PROCESS'] :
                simV2.VisItProcessEngineCommand()
            if command == visit_command['VISIT_COMMAND_SUCCESS'] :
                return True
            if command == visit_command['VISIT_COMMAND_FAILURE'] :
                return False
    return True


def _get_simfilename_from_trace( traceFileName ) :
    """Internal helper function used to extract the name of the simulation file 
       created by libsim from the trace file
       
       Keyword arguments:
       traceFileName : The name of the simulation tracefile
    
    """
    simFile = None
    print "Open tracefile: "+traceFileName
    traceFile = open( traceFileName ,"r")
    print traceFile
    #Find the line containing the sim-filename in the trace file
    for line in traceFile:
        if line.startswith( "Opening sim file " ) :
            entries = line.split(' ') 
            simFile = entries[-1]
            #Remove the endline characters if necessary
            if simFile.endswith("\n") :
                simFile = simFile.rstrip("\n")
    return simFile

