from warp import *
from extpart import *
import visit_engine
import parallel
import numpy
import visit_warpdata



def get_meta_data(cbdata):
    """Callback function used to provide visit with metadata information about the data

       Arguments:
      cbdata : The data object that was handed to VisIt when registering 
                this function as a callback function
    """
    print " get_meta_data: callback called"

    #Allocate a visit metadata object
    metadata = visit_engine.VisIt_SimulationMetaData_alloc()

    #Set the simulation state data
    if cbdata is not None:
        if 'time' in cbdata and 'iteration' in cbdata :
            visit_engine.VisIt_SimulationMetaData_setCycleTime( metadata , cbdata['iteration'] , cbdata['time'] )
        #There are a range of other simulation-specific metadata VisIt_SimulationMetaData_setMode etc.

    #Add the particel data
    if metadata != visit_engine.VISIT_INVALID_HANDLE:
       
       for ex in visit_warpdata.metadata_functions :
           metadata = ex( metadata )

    print " get_meta_data: callback done"
    return metadata


def get_mesh(domain, name, cbdata):
    """Callback function used to send mesh data (e.g., particles) to VisIt.
       The function is used access all meshes the simulation exposes to VisIt.

       VisIt distinguishes between meshes and variables defined on the mesh. 
       We, therefore, need separate callback functions to access the meshes 
       and the variable data defined on the meshes.

       Arguments:
       domain : The domain of the mesh VisIt is requesting
       name : The name of the mesh
       cbdata : The data object that was handed to VisIt when registering 
                this function as a callback function
    """
    print " get_mesh: callback called"

    meshdata = visit_engine.VISIT_INVALID_HANDLE
    for ex in visit_warpdata.meshdata_functions : 
        meshdata = ex( domain, name ,cbdata , meshdata )
        if  meshdata != visit_engine.VISIT_INVALID_HANDLE :
            break
            
    return meshdata


def get_variable(domain, name, cbdata):
    """Callback function used to send variable data (e.g., vx) to VisIt.
    The function is used to access all variables the simulation exposes to VisIt.

    VisIt distinguishes between meshes and variables defined on the mesh. 
    We, therefore, need separate callback functions to access the meshes 
    and the variable data defined on the meshes.

    Arguments:
    domain : The domain of the mesh VisIt is requesting
    name : The name of the mesh
    cbdata : The data object that was handed to VisIt when registering 
             this function as a callback function
    """
    print " get_variable: callback called"
    vardata = visit_engine.VISIT_INVALID_HANDLE
    for ex in visit_warpdata.data_functions :
        vardata = ex( domain, name, cbdata , vardata )
        if vardata != visit_engine.VISIT_INVALID_HANDLE :
            break
    
    print " get_variable: callback done"
    return vardata


def get_domain_list(name, cbdata):
    """Callback function used to get the list of domains handled by this process
    
       In this case we have a single domain per MPI process for all meshes
    
       Arguments:
       name : name of the mesh
       cbdata : The data object that was handed to VisIt when registering 
             this function as a callback function
    """
    
    h = visit_engine.VisIt_DomainList_alloc()
    rank = parallel.get_rank()
    numPE = parallel.number_of_PE()
    if h != VISIT_INVALID_HANDLE:
        hdl = visit_engineVisIt_VariableData_alloc()
        visit_engine.VisIt_VariableData_setDataI(hdl, visit_engine.VISIT_OWNER_VISIT, 1, 1, [rank])
        visit_engine.VisIt_DomainList_setDomains(h, numPE, hdl)
    return h


def control_command_callback(cmd , args , cbdata ) :
    """Callback function used to process user-defined simulation control commands
       issues by the VisIt viewer.
       
       Arguments:
       cmd : The user-defined command issued by the viewer. See control_commands = [...]
       args : Additional arguments provided for the control command.
       cbdata : The data object that was handed to VisIt when registering 
             this function as a callback function
    """
    print " control_command_callback called"
    for ex in visit_warpdata.command_functions :
        status = ex(cmd, args, cbdata )
        if status is True :
            break
    
