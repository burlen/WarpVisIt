from warp import *
from extpart import *
import visit_engine
import parallel
import numpy


######################################################
#        Particle Mesh / Data                        #
######################################################
particle_mesh_name = 'particles'
#Particle mesh variables exposed to VisIt
vx_var_name = 'vx'
vy_var_name = 'vy'
vz_var_name = 'vz'

particle_var_names     = [ 'vx'                                , 'vy'                                 , 'vz' ]
particle_var_types     = [visit_engine.VISIT_VARTYPE_SCALAR    , visit_engine.VISIT_VARTYPE_SCALAR    , visit_engine.VISIT_VARTYPE_SCALAR ]
particle_var_centering = [visit_engine.VISIT_VARCENTERING_NODE , visit_engine.VISIT_VARCENTERING_NODE , visit_engine.VISIT_VARCENTERING_NODE ]

def particle_mesh_metadata( metadata ) :
    """Add the metadata for the particle mesh and associated data to VisIt's metadata object"""
   
    global particle_mesh_name
    global particle_var_names
    global particle_var_types
    global particle_var_centering
    
    # Add particle data
    particlemesh_metadata = visit_engine.VisIt_MeshMetaData_alloc()
    if particlemesh_metadata != visit_engine.VISIT_INVALID_HANDLE:
        # Set the mesh's properties.
        visit_engine.VisIt_MeshMetaData_setName(particlemesh_metadata, particle_mesh_name )
        visit_engine.VisIt_MeshMetaData_setMeshType(particlemesh_metadata, visit_engine.VISIT_MESHTYPE_POINT)
        visit_engine.VisIt_MeshMetaData_setTopologicalDimension(particlemesh_metadata, 0)
        visit_engine.VisIt_MeshMetaData_setSpatialDimension(particlemesh_metadata, 3)
        visit_engine.VisIt_MeshMetaData_setNumDomains(particlemesh_metadata, parallel.number_of_PE())
        visit_engine.VisIt_MeshMetaData_setDomainTitle(particlemesh_metadata, "Domains")
        visit_engine.VisIt_MeshMetaData_setDomainPieceName(particlemesh_metadata, "domain")
        visit_engine.VisIt_MeshMetaData_setNumGroups(particlemesh_metadata, 0)
        visit_engine.VisIt_MeshMetaData_setXUnits(particlemesh_metadata, "m")
        visit_engine.VisIt_MeshMetaData_setYUnits(particlemesh_metadata, "m")
        visit_engine.VisIt_MeshMetaData_setZUnits(particlemesh_metadata, "m")
        visit_engine.VisIt_MeshMetaData_setXLabel(particlemesh_metadata, "x")
        visit_engine.VisIt_MeshMetaData_setYLabel(particlemesh_metadata, "y")
        visit_engine.VisIt_MeshMetaData_setZLabel(particlemesh_metadata, "z")

        visit_engine.VisIt_SimulationMetaData_addMesh(metadata, particlemesh_metadata)

    #Add the particle variable data to the visit metadata object
    for i in xrange(0 , len(particle_var_names) ) :
        
        var_metadata = visit_engine.VisIt_VariableMetaData_alloc()
        if var_metadata != visit_engine.VISIT_INVALID_HANDLE:
            visit_engine.VisIt_VariableMetaData_setName(var_metadata, particle_var_names[i])
            visit_engine.VisIt_VariableMetaData_setMeshName(var_metadata, particle_mesh_name)
            visit_engine.VisIt_VariableMetaData_setType(var_metadata, particle_var_types[i])
            visit_engine.VisIt_VariableMetaData_setCentering(var_metadata, particle_var_centering[i])

            visit_engine.VisIt_SimulationMetaData_addVariable(metadata, var_metadata)

    return metadata


def get_particle_mesh_data( domain, name, cbdata, meshdata):
    """Callback function used to send particle mesh data (e.g., particles) to VisIt.

       Arguments:
       domain : The domain of the mesh VisIt is requesting
       name : The name of the mesh
       cbdata : The data object that was handed to VisIt when registering 
                this function as a callback function
    """
    
    global particle_mesh_name
    #Check for the name of the mesh 
    if name == particle_mesh_name :
        #Allocate the mesh data 
        meshdata = visit_engine.VisIt_PointMesh_alloc()
        if meshdata != visit_engine.VISIT_INVALID_HANDLE :

            xdat = warp.getx()
            ydat = warp.gety()
            zdat = warp.getz()
            rmesh_x = xdat.tolist() #We do not copy the data here, i.e., the simulation should still be the owner of the data not VisIt
            rmesh_y = ydat.tolist()
            rmesh_z = zdat.tolist()
            numParticles = len( rmesh_x )
            hx = visit_engine.VisIt_VariableData_alloc()
            hy = visit_engine.VisIt_VariableData_alloc()
            hz = visit_engine.VisIt_VariableData_alloc()
            visit_engine.VisIt_VariableData_setDataF( hx , visit_engine.VISIT_OWNER_SIM , 1 , numParticles , rmesh_x )
            visit_engine.VisIt_VariableData_setDataF( hy , visit_engine.VISIT_OWNER_SIM , 1 , numParticles , rmesh_y )
            visit_engine.VisIt_VariableData_setDataF( hz , visit_engine.VISIT_OWNER_SIM , 1 , numParticles , rmesh_z )
            visit_engine.VisIt_PointMesh_setCoordsXYZ( meshdata , hx , hy , hz )
            
    return meshdata


def get_particle_mesh_variable(domain, name, cbdata , vardata):
    """Callback function used to send particle mesh variable data (e.g., vx) to VisIt.
    
    Arguments:
    domain : The domain of the mesh VisIt is requesting
    name : The name of the mesh
    vardata : The VisIt metadata object
    """
    global particle_var_names
    if name == particle_var_names[0]  : #vx
        numComponents = 1
        vxData = warp.getvx()
        numTuples = vxData.size
        vardata = visit_engine.VisIt_VariableData_alloc()
        if vardata  != visit_engine.VISIT_INVALID_HANDLE:
            visit_engine.VisIt_VariableData_setDataD(vardata, visit_engine.VISIT_OWNER_SIM, numComponents, numTuples, vxData.tolist())

    elif name == particle_var_names[1] : #vy
        numComponents = 1
        vyData = warp.getvy()
        numTuples = vyData.size
        vardata = visit_engine.VisIt_VariableData_alloc()
        if vardata  != visit_engine.VISIT_INVALID_HANDLE:
            visit_engine.VisIt_VariableData_setDataD(vardata, visit_engine.VISIT_OWNER_SIM, numComponents, numTuples, vyData.tolist())

    elif name == particle_var_names[2] : #Vz
        numComponents = 1
        vzData = warp.getvz()
        numTuples = vzData.size
        vardata = visit_engine.VisIt_VariableData_alloc()
        if vardata  != visit_engine.VISIT_INVALID_HANDLE:
            visit_engine.VisIt_VariableData_setDataD(vardata, visit_engine.VISIT_OWNER_SIM, numComponents, numTuples, vzData.tolist())
            
    return vardata



######################################################
#        Field Mesh / Data                           #
######################################################
field_mesh_name = 'fields'
field_var_names     = [ 'phi'                               , 'rho'                                ]
field_var_types     = [visit_engine.VISIT_VARTYPE_SCALAR    , visit_engine.VISIT_VARTYPE_SCALAR    ]
field_var_centering = [visit_engine.VISIT_VARCENTERING_NODE , visit_engine.VISIT_VARCENTERING_NODE ]
def field_mesh_metadata( metadata ) :
    """Add the metadata for the field mesh and associated data to VisIt's metadata object"""

    global field_mesh_name
    global field_var_names
    global field_var_types
    global field_var_centering
    
    #Add rectilinear field data mesh
    fieldmesh_metadata = visit_engine.VisIt_MeshMetaData_alloc()
    if fieldmesh_metadata != visit_engine.VISIT_INVALID_HANDLE:
        # Set the mesh's properties.
        visit_engine.VisIt_MeshMetaData_setName(fieldmesh_metadata, field_mesh_name )
        visit_engine.VisIt_MeshMetaData_setMeshType(fieldmesh_metadata, visit_engine.VISIT_MESHTYPE_RECTILINEAR)
        visit_engine.VisIt_MeshMetaData_setTopologicalDimension(fieldmesh_metadata, 3)
        visit_engine.VisIt_MeshMetaData_setSpatialDimension(fieldmesh_metadata, 3)
        visit_engine.VisIt_MeshMetaData_setNumDomains(fieldmesh_metadata, parallel.number_of_PE())
        visit_engine.VisIt_MeshMetaData_setDomainTitle(fieldmesh_metadata, "Domains")
        visit_engine.VisIt_MeshMetaData_setDomainPieceName(fieldmesh_metadata, "domain")
        visit_engine.VisIt_MeshMetaData_setNumGroups(fieldmesh_metadata, 0)
        visit_engine.VisIt_MeshMetaData_setXUnits(fieldmesh_metadata, "m")
        visit_engine.VisIt_MeshMetaData_setYUnits(fieldmesh_metadata, "m")
        visit_engine.VisIt_MeshMetaData_setZUnits(fieldmesh_metadata, "m")
        visit_engine.VisIt_MeshMetaData_setXLabel(fieldmesh_metadata, "x")
        visit_engine.VisIt_MeshMetaData_setYLabel(fieldmesh_metadata, "y")
        visit_engine.VisIt_MeshMetaData_setZLabel(fieldmesh_metadata, "z")

        visit_engine.VisIt_SimulationMetaData_addMesh(metadata, fieldmesh_metadata)
        
    #Add the field variable data to the visit metadata object
    for i in xrange(0 , len(field_var_names) ) :
        
        var_metadata = visit_engine.VisIt_VariableMetaData_alloc()
        if var_metadata != visit_engine.VISIT_INVALID_HANDLE:
            visit_engine.VisIt_VariableMetaData_setName(var_metadata, field_var_names[i])
            visit_engine.VisIt_VariableMetaData_setMeshName(var_metadata, field_mesh_name)
            visit_engine.VisIt_VariableMetaData_setType(var_metadata, field_var_types[i])
            visit_engine.VisIt_VariableMetaData_setCentering(var_metadata, field_var_centering[i])

            visit_engine.VisIt_SimulationMetaData_addVariable(metadata, var_metadata)

    return metadata


def get_field_mesh_data( domain, name, cbdata , meshdata):
    """Callback function used to send field  mesh data (e.g., particles) to VisIt.

       Arguments:
       domain : The domain of the mesh VisIt is requesting
       name : The name of the mesh
       cbdata : The data object that was handed to VisIt when registering 
                this function as a callback function
    """
    
    global field_mesh_name
    if name == field_mesh_name :
        #Allocate the mesh data 
        meshdata = visit_engine.VisIt_RectilinearMesh_alloc()
        if meshdata != visit_engine.VISIT_INVALID_HANDLE :

            dx = w3d.dx
            dy = w3d.dy
            dz = w3d.dz
            rmesh_x = numpy.arange( start=w3d.xmmin , stop= (w3d.xmmax + (dx/2.0)) , step=dx ).tolist() #The upper bound is modified to ensure that arange includes the upper bound
            rmesh_y = numpy.arange( start=w3d.ymmin , stop= (w3d.ymmax + (dy/2.0)) , step=dy ).tolist()
            rmesh_z = numpy.arange( start=w3d.zmmin , stop= (w3d.zmmax + (dz/2.0)) , step=dz ).tolist()
            nx = len( rmesh_x ) #Same as w3d.nx
            ny = len( rmesh_y ) #Same as w3d.ny
            nz = len( rmesh_z ) #Same as w3d.nz
            hx = visit_engine.VisIt_VariableData_alloc()
            hy = visit_engine.VisIt_VariableData_alloc()
            hz = visit_engine.VisIt_VariableData_alloc()
            visit_engine.VisIt_VariableData_setDataF( hx , visit_engine.VISIT_OWNER_SIM , 1 , nx , rmesh_x )
            visit_engine.VisIt_VariableData_setDataF( hy , visit_engine.VISIT_OWNER_SIM , 1 , ny , rmesh_y )
            visit_engine.VisIt_VariableData_setDataF( hz , visit_engine.VISIT_OWNER_SIM , 1 , nz , rmesh_z )
            visit_engine.VisIt_RectilinearMesh_setCoordsXYZ( meshdata , hx , hy , hz )
            
    return meshdata
    
    
def get_field_mesh_variable(domain, name, cbdata, vardata):
    """Callback function used to send field mesh variable data (e.g., vx) to VisIt.
    
    Arguments:
    domain : The domain of the mesh VisIt is requesting
    name : The name of the mesh
    vardata : The VisIt metadata object
    """
    
    global field_var_names
    if name == field_var_names[0] : #phi
        numComponents = 1
        numTuples = warp.getphi().size
        phiData = numpy.reshape( warp.getphi() , numTuples , order='F' ).tolist() #NOTE: We need to ensure fortran ordering
        vardata = visit_engine.VisIt_VariableData_alloc()
        if vardata  != visit_engine.VISIT_INVALID_HANDLE:
            visit_engine.VisIt_VariableData_setDataD(vardata, visit_engine.VISIT_OWNER_COPY, numComponents, numTuples, phiData)
        else :
            print "Allocation for variable data object failed"

    elif name == field_var_names[1] : #rho
        numComponents = 1
        numTuples = warp.getrho().size
        rhoData = numpy.reshape( warp.getrho() , numTuples , order='F' ).tolist() #NOTE: We need to ensure fortran ordering
        vardata = visit_engine.VisIt_VariableData_alloc()
        if vardata  != visit_engine.VISIT_INVALID_HANDLE:
            visit_engine.VisIt_VariableData_setDataD(vardata, visit_engine.VISIT_OWNER_SIM, numComponents, numTuples, rhoData )
            
    return vardata


######################################################
#        CSG Mesh /Data                              #
######################################################
csg_mesh_name = 'csg'
#Geometry settings
piperad_outer = 3.445e-2 + 0.3
piperad_inner = piperad_outer * 0.8 
pipelength = 2.0 #4.

def csg_mesh_metadata( metadata ) :
    """Add the metadata for the csg mesh and associated data to VisIt's metadata object"""
    
    #Add rectilinear field data mesh
    csgmesh_metadata = visit_engine.VisIt_MeshMetaData_alloc()
    if csgmesh_metadata != visit_engine.VISIT_INVALID_HANDLE:
        # Set the mesh's properties.
        visit_engine.VisIt_MeshMetaData_setName(csgmesh_metadata, csg_mesh_name )
        visit_engine.VisIt_MeshMetaData_setMeshType(csgmesh_metadata, visit_engine.VISIT_MESHTYPE_CSG)
        visit_engine.VisIt_MeshMetaData_setTopologicalDimension(csgmesh_metadata, 3)
        visit_engine.VisIt_MeshMetaData_setSpatialDimension(csgmesh_metadata, 3)
        visit_engine.VisIt_MeshMetaData_setNumDomains(csgmesh_metadata, parallel.number_of_PE())
        visit_engine.VisIt_MeshMetaData_setDomainTitle(csgmesh_metadata, "Domains")
        visit_engine.VisIt_MeshMetaData_setDomainPieceName(csgmesh_metadata, "domain")
        visit_engine.VisIt_MeshMetaData_setNumGroups(csgmesh_metadata, 0)
        visit_engine.VisIt_MeshMetaData_setXUnits(csgmesh_metadata, "m")
        visit_engine.VisIt_MeshMetaData_setYUnits(csgmesh_metadata, "m")
        visit_engine.VisIt_MeshMetaData_setZUnits(csgmesh_metadata, "m")
        visit_engine.VisIt_MeshMetaData_setXLabel(csgmesh_metadata, "x")
        visit_engine.VisIt_MeshMetaData_setYLabel(csgmesh_metadata, "y")
        visit_engine.VisIt_MeshMetaData_setZLabel(csgmesh_metadata, "z")

        visit_engine.VisIt_SimulationMetaData_addMesh(metadata, csgmesh_metadata)
    #We do not have any direct data on the CSG mesh right now
    
    return metadata
    
    
def get_field_mesh_data( domain, name, cbdata , meshdata):
    """Callback function used to send csg mesh data (e.g., particles) to VisIt.

       Arguments:
       domain : The domain of the mesh VisIt is requesting
       name : The name of the mesh
       cbdata : The data object that was handed to VisIt when registering 
                this function as a callback function
    """
    
    global csg_mesh_name
    global piperad_inner
    global piperad_outer
    global pipelength
    
    if name == csg_mesh_name :
        
        global piperad
        global pipelength
        zmin = top.zbeam - pipelength/2.
        zmax = top.zbeam + pipelength/2.
        csg_extends_min = [-piperad_outer , -piperad_outer , zmin]
        csg_extends_max = [ piperad_outer ,  piperad_outer , zmax]
        csg_bound_types = [ visit_engine.VISIT_CSG_CYLINDER_PNLR , \
                            visit_engine.VISIT_CSG_CYLINDER_PNLR ]
        csg_bound_coeffs = [ 0. , 0. , zmin , 0. , 0. , 1. , pipelength , piperad_outer , \
                             0. , 0. , zmax , 0. , 0. , 1. , pipelength , piperad_inner ]
        csg_region_operators = [ visit_engine.VISIT_CSG_INNER , visit_engine.VISIT_CSG_OUTER , visit_engine.VISIT_CSG_INTERSECT ]
        csg_leftids =  [ 0 ,  1, 0]
        csg_rightids = [-1 , -1, 1]
        csg_zonelist = [2]
       
        #Allocate the mesh data 
        meshdata = visit_engine.VisIt_CSGMesh_alloc()
        if meshdata != visit_engine.VISIT_INVALID_HANDLE :
            
            cbt = visit_engine.VisIt_VariableData_alloc()
            visit_engine.VisIt_VariableData_setDataI( cbt , visit_engine.VISIT_OWNER_COPY , 1 , len(csg_bound_types) , csg_bound_types)
            visit_engine.VisIt_CSGMesh_setBoundaryTypes( meshdata , cbt )
            
            cbc = visit_engine.VisIt_VariableData_alloc()
            visit_engine.VisIt_VariableData_setDataF( cbc , visit_engine.VISIT_OWNER_COPY , 1 , len(csg_bound_coeffs) , csg_bound_coeffs)
            visit_engine.VisIt_CSGMesh_setBoundaryCoeffs( meshdata , cbc )
            
            visit_engine.VisIt_CSGMesh_setExtents( meshdata , csg_extends_min , csg_extends_max )
            
            cro = visit_engine.VisIt_VariableData_alloc()
            visit_engine.VisIt_VariableData_setDataI( cro , visit_engine.VISIT_OWNER_COPY , 1 , len(csg_region_operators) , csg_region_operators)
            cli = visit_engine.VisIt_VariableData_alloc()
            visit_engine.VisIt_VariableData_setDataI( cli , visit_engine.VISIT_OWNER_COPY , 1 , len(csg_leftids) , csg_leftids)
            cri = visit_engine.VisIt_VariableData_alloc()
            visit_engine.VisIt_VariableData_setDataI( cri , visit_engine.VISIT_OWNER_COPY , 1 , len(csg_rightids) , csg_rightids)
            visit_engine.VisIt_CSGMesh_setRegions( meshdata , cro , cli, cri )
            
            czl = visit_engine.VisIt_VariableData_alloc()
            visit_engine.VisIt_VariableData_setDataI( czl , visit_engine.VISIT_OWNER_COPY , 1 , len(csg_zonelist) , csg_zonelist)
            visit_engine.VisIt_CSGMesh_setZonelist( meshdata , czl )

    return meshdata


######################################################
#        Expressions/Derived Quantities              #
######################################################
expression_names = [ 'vvec'        , 'energy'                      , 'const' ]
expression_defs  = [ "{vx, vy, vz}", "sqrt(vx*vx + vy*vy + vz*vz)" , "coord(csg)[0]*0" ]
expression_types = [ visit_engine.VISIT_VARTYPE_VECTOR, visit_engine.VISIT_VARTYPE_SCALAR, visit_engine.VISIT_VARTYPE_SCALAR ]

def expression_metadata( metadata ) :
    """Add all expressions to visit's metadata object"""

    global expression_names
    global expression_defs
    global expression_types

    for ei in xrange(0 , len(expression_names) ) :
        temp = visit_engine.VisIt_ExpressionMetaData_alloc()
        if temp != visit_engine.VISIT_INVALID_HANDLE:
            visit_engine.VisIt_ExpressionMetaData_setName(temp, expression_names[ei] )
            visit_engine.VisIt_ExpressionMetaData_setDefinition(temp, expression_defs[ei] )
            visit_engine.VisIt_ExpressionMetaData_setType(temp, expression_types[ei] )

            visit_engine.VisIt_SimulationMetaData_addExpression(metadata, temp)

    return metadata



######################################################
#        Control Commands                            #
######################################################
#List of commands the viewer may use
control_commands = ["halt" , "step" , "run", "run_without_update" , "end"]
def control_commands_metadata( metadata ) :
    """Add all control commands to visit's metadata object"""
    for c in control_commands:
        cmd = visit_engine.VisIt_CommandMetaData_alloc()
        if cmd != visit_engine.VISIT_INVALID_HANDLE:
            visit_engine.VisIt_CommandMetaData_setName(cmd, c)
            visit_engine.VisIt_SimulationMetaData_addGenericCommand(metadata, cmd)

    return metadata

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
    status = True
    if cmd == "halt" :
        cbdata['run_mode'] = 0
        visit_engine.VisItUpdatePlots()
        print "Halted the simulation"
    elif cmd == "step":
        for i in xrange( 0 , cbdata['step_size'] ) :
            print "-------------------------------"+str(i)
            step() 
            cbdata ['iteration'] = top.it
            cbdata ['time'] = top.time
            visit_engine.VisItTimeStepChanged()
        visit_engine.VisItUpdatePlots()
        print "Ran "+str(cbdata['step_size'])+" simulation step(s)"
    elif cmd == "run" : #run
        cbdata['run_mode'] = 1
        print "Run the simulation"
    elif cmd == "run_without_update" : #run_without_update
        cbdata['run_mode'] = 2
        print "Run the simulation without updating the visualization"
    elif cmd == "end" :
        cbdata['run_mode'] = 3
        print "End the simulation main loop"
    else :
        status = False
    return status
    


######################################################
# Define the functions that need to be exposed to    #
# visit_datainterface to define the data that VisIt  #
# should know about.                                 #
######################################################
#List of functions defining VisIt metadata
metadata_functions = [ particle_mesh_metadata , field_mesh_metadata , csg_mesh_metadata  , expression_metadata , control_commands_metadata ]

#List of mesh data functions
meshdata_functions = [ get_particle_mesh_data , get_field_mesh_data, get_particle_mesh_data ]

#Data callback functions
data_functions = [ get_particle_mesh_variable , get_field_mesh_variable ]

#Command callback functions
command_functions =  [ control_command_callback ]