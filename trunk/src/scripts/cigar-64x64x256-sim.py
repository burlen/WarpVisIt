from warp import *
from extpart import *

def UpdatePlots():
    """
    Return true if the current time step should be visualized.
    """
    # this saves an image every ten steps
    return (int(warp.top.it) % int(10)) == int(0)

def Finished():
    """
    Return true when the simulation should no longer run
    """
    return warp.top.it >= 100

def Initialize():
    """
    Example of a how an in-situ run is configured.
    This funcion must be named Configure. It is
    passed to the in-situ driver via environment
    variable name WarpConfigFile set with the path
    to this file.
    """
    # --- Set four-character run id, comment lines, user's name.
    top.pline2   = "Example 3D beam in a FODO lattice"
    top.pline1   = "S-G cigar beam. 64x64x256"
    top.runmaker = "David P. Grote"

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
    addnewquad(zs= -quadlen/2.,
           ze= +quadlen/2.,
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
    return
