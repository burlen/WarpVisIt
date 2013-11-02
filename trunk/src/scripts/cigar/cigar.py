from warp import *
from extpart import *
import os

#-----------------------------------------------------------------------------
def LoadRenderScripts(scriptRoot):
    """
    This function returns a dictionary of rendering scripts
    whose key is a descriptive string. This dictionary will
    be used when our Render function returns a list of scripts
    to run. the scriptRoot argument gives the path where scripts
    are stored.
    """
    # user supplied rendering scripts
    # script names (keys) are used by the Render funciton
    # to select the desired script
    renderingScripts = {
        'max(v(x,y))' : 'render-max-v.py',
        'binning bv'  : 'render-binning-bv.py',
        'scatter bv'  : 'render-scatter-bv.py',
        'volume phi'  : 'render-volume-phi.py',
        'particle v'  : 'render-particle-v.py',
        }

    for key,fileName in renderingScripts.iteritems():
        f = open(os.path.join(scriptRoot,fileName))
        code = f.read()
        f.close()
        renderingScripts[key] = code

    return renderingScripts


#-----------------------------------------------------------------------------
def GetActiveRenderScripts():
    """
    If the current time step should be visualized return a list
    of keys naming which rendering scripts should run. If the
    list is empty then nothing will be rendered this time step.
    The script dictionary is created by LoadRenderScripts.
    """
    scripts = []

    # some very contrived examples
    # this just shows the flexibility
    # the point is rendering can be triggered for
    # any plot under any condition

    # after M iterations every N iterations plot...
    #if (warp.top.it >= 10) and ((warp.top.it % 8) == 0):
    #    scripts.append('particle v')
    #    scripts.append('volume phi')

    # every N interations plot ...
    #if ((warp.top.it % 5) == 0):
    #    scripts.append('max(v(x,y))')
    #    scripts.append('num b vs. v')

    # always plot
    scripts.append('scatter bv')
    scripts.append('binning bv')
    scripts.append('particle v')
    scripts.append('volume phi')
    scripts.append('max(v(x,y))')

    return scripts

#-----------------------------------------------------------------------------
def Advance():
    """Advance the simulation one time step."""
    warp.step()

#-----------------------------------------------------------------------------
def Continue():
    """
    Return false when the simulation should no longer run.
    """
    # adjust this to take as many steps as you need
    return warp.top.it <= 4


#-----------------------------------------------------------------------------
def Finalize():
    """
    shutdown the simulation, cleanup etc.
    """
    pass

#-----------------------------------------------------------------------------
def Initialize():
    """
    Setup IC and start the simulation, but don't run it yet.
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
