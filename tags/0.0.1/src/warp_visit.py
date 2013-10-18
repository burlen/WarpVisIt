from warp import *
from extpart import *
import parallel

import visit_engine
import visit_viewer
import visit_datainterface
import visit_warpdata


"""Additional data passed to the VisIt callback functions. This is used to share data between the viewer callbacks and the simulation.
   iteration : Simulation iteration index    (updated by the simulation)
   time : Simulation time   (updated by the simulation)
   run_mode : 0=Wait for VisIt viewer control input.
              1=Run the simualtion and update visualization. I.e., timeout if no input from VisIt viewer is found.
              2=Run the simualtion without updating the visualization. I.e., timeout if no input from VisIt viewer is found.
              3=Terminate the main simulation loop
   step_size : Number of simulations steps to be executed when stepping throuhgh the simulation
   nstep : Maximum number of simulation steps to be executed. Set to -1 to indicate no upper limit.
   sim_file : The libsim simulation file which allows VisIt to connect to the simulation. (automatically set in main)
"""
callback_userdata = {'iteration':0, 'time':0, 'run_mode':0, 'step_size':2, 'nstep':1000 , 'sim_file':None}


def main(argv=None):

    global callback_userdata

    #Brief description of the simulation
    sim_title = "Example 3D beam in a FODO lattice"
    sim_comment = "S-G cigar beam. 64x64x256"
    sim_runmaker = "David P. Grote"

    #Settings for the VisIt viewer
    visit_viewer_rank = 0             #Which MPI rank should VisIt's viewer run on
    visit_viewer_rank_shared = True  #Does the VisIt viewer share the process (i.e., we need to use multiprocessing) or does the simulation run in a sepearte process
    visit_viewer_show_window = True   #Should a viewer window be launched
    visit_viewer_debug = -1           #Integer indicating the debug level for the VisIt viewer. -1 (disable), 1,2,3,4,5 (enable debug level)
    visit_viewer_args = []            #Optional launch arguments for the VisIt viewer
    visit_viewer_external = False     #Use an external VisIt viewer rather than one that we start ourselfs

    #Initalize the simulation
    setup_simulation(sim_title , sim_comment, sim_runmaker)
    print "Setup simulation complete"

    #Initalize the visit comute engine
    #Additional input parameters for visit_engine.initalize(...) left as default here are:
    #sim_path = os.getcwd() #Path were the simulation is running.
    #visitTraceFile = None  #Define the file were the output trace of libsim should be stored.
    print "Initalize the VisIt libsim engine"
    simFile = visit_engine.initalize( sim_comment = (sim_runmaker+": "+sim_title+"; "+sim_comment) , visitTraceFile="/global/u2/o/oruebel/Warp_Test/w_simTrace.txt" ) 
    print "Initalized VisIt libsim engine"
    
     #Execute the simulation main function
    callback_userdata['run_mode'] = 0
    #import multiprocessing 
    #tsim = multiprocessing.Process(target=sim_main)
    #tsim.daemon = True #By making the viewer process a daemon allows us to call exit
    #tsim.start()
    
    #Wait for a visit viewer to connect
    if simFile is not None :
        print simFile
        callback_userdata['sim_file'] = simFile
        #Start the VisIt viewer and establish the initial connection
        if not visit_viewer_external :
            print "Starting VisIt Viewer"
            visit_viewer.start_visit_viewer( simFile=simFile , vis_function=run_visualization, \
                              visit_viewer_rank=visit_viewer_rank , visit_viewer_rank_shared=visit_viewer_rank_shared , \
                              visit_viewer_show_window=visit_viewer_show_window , visit_viewer_debug=visit_viewer_debug, \
                              visit_viewer_args=visit_viewer_args )
#        print "VisIt engine waiting for the viewer to connect"
#        if True : # not visit_viewer_external or (visit_wait_for_external_viewer and visit_viewer_external) :
#            visit_viewer_connected = visit_engine.waitForViewer(visit_datainterface, callback_userdata)
#            if visit_viewer_connected :
#                print "VisIt Engine: connection complete"
#            else:
#                print "VisIt Engine: connection with viewer failed"
#        else :
#             callback_userdata['run_mode'] = 1
    else:
        print "Initalization of VisIt simV2 failed"
        exit()

    sim_main()

    print "Simulation done"


def sim_main() :
    """Main loop for the simulation which reacts to input by the VisIt viewer"""

    global callback_userdata
    nStep=callback_userdata['nstep'] 
    step_size=callback_userdata['step_size']
    error = False  #Did we encounter an error
    status = 0     #VisIt input status

    #Run the simulation as long as:
    # 1) We have not reached the maximum number of simulation steps
    # 2) The VisIt viewer has not told simulation explicitly to terminate the main loop
    while (nStep<0 or callback_userdata['iteration']<nStep) and callback_userdata['run_mode'] != 3 :

        #Check for input controls by the VisIt viewer
        if parallel.get_rank() == 0 :
            if callback_userdata['run_mode']== 0 :
                status = visit_engine.VisItDetectInput(True , -1)
            else :
                status = visit_engine.VisItDetectInputWithTimeout(False , 1 , -1 )
        visit_engine.broadcast_int( status , 0)
        #print "Status: "+str(status)

        #Depending on the response by the VisIt viewer take one of the following actions
        #An error occured 
        if status <= -1 : 
            print "Error"
            error=True
            break

        #Simulate a single timestep. (Ok - Timed out)
        elif status == 0 : 
            print "Simulate "+str(step_size)+"timestep(s)"
            for i in xrange( 0 , step_size ) :
                step() 
                callback_userdata['iteration'] = top.it
                callback_userdata['time'] = top.time
                visit_engine.VisItTimeStepChanged() #Notify VisIt that the timestep has changed
            #Notify VisIt that the simulation data has changed
            if callback_userdata['run_mode'] != 2 :
                visit_engine.VisItUpdatePlots()

        #Establish connection to another VisIt viewer. (Listen socket input)
        elif status == 1:
            print "Another viewer tries to connect"
            if visit_engine.VisItAttemptToCompleteConnection() :
                visit_engine.connect_callbacks(visit_datainterface, callback_userdata)

        #Process a command for the VisIt compute engine (Engine socket input)
        elif status == 2:
            exStat = visit_engine.process_visit_command()
            if not exStat :
                print "Closing connection to VisIt"
                visit_engine.VisItDisconnect()
            else :
                pass
                #print "Processed VisIt engine command"

        #Console socket input
        elif status == 3:
            raise NotImplementedError("Console socket input is currently not handled by the simulation")

    #Finalize the visualitation (close the trace file etc.)
    visit_engine.finalize()
    
    #Return the error status
    return error



def run_visualization(simFile) :

    print dir(visit_viewer)
    AnnotationAtts = visit_viewer.visit.GetAnnotationAttributes()
    AnnotationAtts.axes3D.visible = 0
    AnnotationAtts.axes3D.bboxFlag = 0
    AnnotationAtts.userInfoFlag = 0
    visit_viewer.visit.SetAnnotationAttributes(AnnotationAtts)
    
    RenderingAtts = visit_viewer.visit.GetRenderingAttributes()
    RenderingAtts.antialiasing = 1  
    RenderingAtts.scalableActivationMode = RenderingAtts.Always  # Never, Always, Auto  
    RenderingAtts.doShadowing = 1
    RenderingAtts.shadowStrength = 0.5
    visit_viewer.visit.SetRenderingAttributes(RenderingAtts)
    
    light0 = visit_viewer.visit.LightAttributes()
    light0.enabledFlag = 1
    light0.type = light0.Camera  # Ambient, Object, Camera
    #light0.direction = (-0.363, -0.508, -0.781)
    light0.direction = ( 0.2, -0.2, -0.8)
    light0.color = (255, 255, 255, 255)
    light0.brightness = 1
    visit_viewer.visit.SetLight(0, light0)

    import time
    #Add plots 
    #visit.AddPlot( "Mesh", visit_datainterface.particle_mesh_name )
    visit_viewer.visit.AddPlot( "Pseudocolor", "vx" )
    #visit.AddPlot( "Pseudocolor", visit_datainterface.expression_names[1] )
    #visit.AddPlot( "Vector" , visit_datainterface.vvec_exp_name )
    #visit.AddPlot( "Mesh", visit_datainterface.field_mesh_name )
    #AddPlot( "Pseudocolor", visit_datainterface.rho_var_name )
    PseudocolorAtts = visit_viewer.visit.PseudocolorAttributes()
    PseudocolorAtts.legendFlag = 1
    PseudocolorAtts.lightingFlag = 1
    PseudocolorAtts.pointSize = 0.005
    PseudocolorAtts.pointType = PseudocolorAtts.Icosahedron  # Box, Axis, Icosahedron, Point, Sphere
    PseudocolorAtts.opacity = 1.0
    PseudocolorAtts.colorTableName = "hot_desaturated"
    PseudocolorAtts.pointSizePixels = 2
    PseudocolorAtts.opacityType = PseudocolorAtts.Explicit  # Explicit, ColorTable
    visit_viewer.visit.SetPlotOptions(PseudocolorAtts)
    
    
    visit_viewer.visit.AddPlot( "Pseudocolor", visit_warpdata.expression_names[2] )
    visit_viewer.visit.AddOperator("Clip", 0)
    ClipAtts = visit_viewer.visit.ClipAttributes()
    ClipAtts.funcType = ClipAtts.Plane  # Plane, Sphere
    ClipAtts.plane1Status = 1
    ClipAtts.plane2Status = 0
    ClipAtts.plane3Status = 0
    ClipAtts.plane1Origin = (-0.15, 0, 0)
    ClipAtts.plane1Normal = (1, 0, 0)
    ClipAtts.planeInverse = 1
    ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1  # None, Plane1, Plane2, Plane3
    ClipAtts.center = (0, 0, 0)
    ClipAtts.radius = 1
    ClipAtts.sphereInverse = 0
    visit_viewer.visit.SetOperatorOptions(ClipAtts, 0)
    
    print "VisIt Viewer: added plot"
    visit_viewer.visit.DrawPlots()
    time.sleep(20)
    #SaveWindow()
    print "VisIt Viewer: plot drawn"#time.sleep(3)
    #time.sleep(3)
    for i in xrange( 0 , 10 ) :
        visit_viewer.visit.SendSimulationCommand('localhost', simFile, 'step')
        visit_viewer.visit.DrawPlots()
        #SaveWindow()
    visit_viewer.visit.SendSimulationCommand('localhost', simFile, 'run_without_update')
    time.sleep(20)
    visit_viewer.visit.SendSimulationCommand('localhost', simFile, 'halt')
    visit_viewer.visit.DrawPlots()
    #SaveWindow()
    time.sleep(100)
    #SendSimulationCommand('localhost', simFile, 'halt')
    #DrawPlots()
    #SaveWindow()
    #time.sleep(60)
    visit_viewer.visit.SendSimulationCommand('localhost', simFile, 'end')
    exit(0)



def setup_simulation(sim_title="", sim_comment="", sim_runmaker=""  ) :
    """Setup the simulation 
    
       Optional Arguments:
       sim_title : Title of the simulation
       sim_comment : Comment for the simulation
       sim_runmaker : The user running the simulation
    """
    
    # --- Set four-character run id, comment lines, user's name.
    top.pline2   = sim_title
    top.pline1   = sim_comment
    top.runmaker = sim_runmaker

    # --- Invoke setup routine - it is needed to created a cgm file for plots
    setup()

    # --- Create the beam species
    beam = Species(type=Potassium,charge_state=+1,name="Beam species")

    # --- Set input parameters describing the beam, 72 to 17.
    beam.b0       = 15.358933450767e-3
    beam.a0       =  8.6379155933081e-3
    beam.x0       = 3.*mm
    beam.emit     = 51.700897052724e-6
    beam.ap0      = 0.e0
    beam.bp0      = 0.e0
    beam.ibeam    = 2.e-03
    beam.vbeam    = 0.e0
    beam.ekin     = 80.e3
    beam.aion     = beam.type.A
    beam.zion     = beam.charge_state
    top.lrelativ = false
    top.derivqty()
    beam.vthz     = .5e0*beam.vbeam*beam.emit/sqrt(beam.a0*beam.b0) # Vthz ~ Vthperp

    # +++ Set up arrays describing lattice.
    # --- Set temp variables.
    hlp     = 36.0e-2  # half lattice period length
    piperad = 3.445e-2 # pipe radius
    quadlen = 11.e-2   # quadrupole length
    gaplen = 4.*cm
    rodlen = quadlen + gaplen
    dbdx    = .949/quadlen

    # --- Set general lattice variables.
    top.tunelen   = 2.e0*hlp
    env.zl        = -hlp*2
    env.zu        = -env.zl
    env.dzenv     = top.tunelen/100.e0

    # --- Set up quadrupoles
    addnewquad(zs=    - quadlen/2.,
           ze=    + quadlen/2.,
           db=-dbdx,ap=piperad)
    addnewquad(zs=hlp - quadlen/2.,
           ze=hlp + quadlen/2.,
           db=+dbdx,ap=piperad)
    addnewquad(zs=2.*hlp - quadlen/2.,
           ze=2.*hlp + quadlen/2.,
           db=-dbdx,ap=piperad)
    top.zlatstrt  = 0.
    top.zlatperi  = 2.e0*hlp


    # +++ Set input parameters describing the 3d simulation.
    w3d.nx = 64/2
    w3d.ny = 64/2
    w3d.nz = 256/2
    steps_p_perd = 50
    top.dt = (top.tunelen/steps_p_perd)/beam.vbeam

    # --- Set to finite beam.
    top.pbound0  = top.pboundnz = periodic
    top.pboundxy = absorb
    w3d.xmmin = -piperad
    w3d.xmmax =  piperad
    w3d.ymmin = -piperad
    w3d.ymmax =  piperad
    w3d.zmmin = -hlp*2
    w3d.zmmax = +hlp*2
    top.prwall = piperad

    # --- Set pulse length.
    beam.zimin = w3d.zmmin*.95/2.
    beam.zimax = w3d.zmmax*.95/2.

    # --- Load Semi-Gaussian cigar beam.
    top.npmax = 20000
    w3d.distrbtn = "semigaus"
    w3d.cigarld = true
    w3d.xrandom = "digitrev"
    w3d.vtrandom = "digitrev"
    w3d.vzrandom = "digitrev"
    w3d.ldprfile = "polar"
    w3d.cylinder = false
    top.straight = .8

    # --- set up field solver
    w3d.l4symtry = true
    w3d.bound0 = periodic
    w3d.boundnz = periodic
    w3d.boundxy = dirichlet

    solver = MultiGrid3D()
    registersolver(solver)

    pipe = ZCylinderOut(piperad,4.,voltage=0.)
    installconductors(pipe,dfill=largepos)

    # --- Run the envelope solver to provide data used to initialize particles.
    package("env")
    generate()
    step()

    # --- Generate the PIC code (allocate storage, load ptcls, t=0 plots, etc.).
    package("w3d")
    generate()

#def visualize(lwrite=0):
   # pipe.zcent = top.zbeam
   # pipe.createdxobject(rend=piperad+0.3*cm,phimin=-pi,phimax=0.,
   #                     color=[0.4,0.4,0.4],fillinends=1,normalsign=-1,close=1)
   # dd = Opyndx.DXCollection(pipe)
   #
   # iqmin = int((pipe.zcent - pipe.length/2.)/hlp)
   # iqmax = int((pipe.zcent + pipe.length/2.)/hlp)
   # for iq in range(iqmin,iqmax + 1):
   #   quad = ZCylinderOut(piperad-1.*um,quadlen,zcent=iq*hlp)
   #   listofallconductors.remove(quad)
   #   quad.createdxobject(rend=piperad,phimin=-pi,phimax=0.,
   #                       color=[0.6,0.4,0.4],fillinends=1,normalsign=-1,close=0)
   #   dd.addobject(quad)
   #
   # vx = getvx()
   # dotsize = 0.6e-3
   # ipstep = 2
   # ppp,colormap = Opyndx.viewparticles(
   #                       getx()[::ipstep],gety()[::ipstep],getz()[::ipstep],
   #                       vx[::ipstep]/18000.,color='auto',colorbar=1,
   #                       ratio=1.,
   #                       type=0.2,
   #                       scale=[1./dotsize,1./dotsize,1./dotsize],
   #                       display=0)
   # dd.addobject(ppp)
   #
   # camera = Opyndx.DXCamera([0.,0.,top.zbeam],
   #                          [.005,0.06,-1.0+top.zbeam],
   #                          up=[1,0,0],width=0.2)
   #
   # ambient = Opyndx.DXAmbientLight([.6,.6,.6])
   # light = Opyndx.DXLight([0.0,-5.,0.],[1.,1.,1.])
   # dd.addobject(ambient)
   # dd.addobject(light)
   #
   # if lwrite:
   #   Opyndx.DXWriteImage('image%04d.tiff'%(top.it),dd,camera=camera,
   #                       format='tiff')
   # else:
   #   Opyndx.DXImage(dd,camera)

#visualize()


#Initalize the visit viewer used to control the visualization
#import visit_viewer


#def makemovie(nsteps=50,isteps=2):
#  if top.it == 0: visualize(lwrite=1)
#  while top.it < nsteps:
#    step(isteps)
#    visualize(lwrite=1)



if __name__ == "__main__":
    main()

