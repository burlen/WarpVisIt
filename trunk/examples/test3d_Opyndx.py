from warp import *
from extpart import *
import Opyndx

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

def visualize(lwrite=0):
    pipe.zcent = top.zbeam
    pipe.createdxobject(rend=piperad+0.3*cm,phimin=-pi,phimax=0.,
                        color=[0.4,0.4,0.4],fillinends=1,normalsign=-1,close=1)
    dd = Opyndx.DXCollection(pipe)

    iqmin = int((pipe.zcent - pipe.length/2.)/hlp)
    iqmax = int((pipe.zcent + pipe.length/2.)/hlp)
    for iq in range(iqmin,iqmax + 1):
      quad = ZCylinderOut(piperad-1.*um,quadlen,zcent=iq*hlp)
      listofallconductors.remove(quad)
      quad.createdxobject(rend=piperad,phimin=-pi,phimax=0.,
                          color=[0.6,0.4,0.4],fillinends=1,normalsign=-1,close=0)
      dd.addobject(quad)

    vx = getvx()
    dotsize = 0.6e-3
    ipstep = 2
    ppp,colormap = Opyndx.viewparticles(
                          getx()[::ipstep],gety()[::ipstep],getz()[::ipstep],
                          vx[::ipstep]/18000.,color='auto',colorbar=1,
                          ratio=1.,
                          type=0.2,
                          scale=[1./dotsize,1./dotsize,1./dotsize],
                          display=0)
    dd.addobject(ppp)

    camera = Opyndx.DXCamera([0.,0.,top.zbeam],
                             [.005,0.06,-1.0+top.zbeam],
                             up=[1,0,0],width=0.2)

    ambient = Opyndx.DXAmbientLight([.6,.6,.6])
    light = Opyndx.DXLight([0.0,-5.,0.],[1.,1.,1.])
    dd.addobject(ambient)
    dd.addobject(light)

    if lwrite:
      Opyndx.DXWriteImage('image%04d.tiff'%(top.it),dd,camera=camera,
                          format='tiff')
    else:
      Opyndx.DXImage(dd,camera)

visualize(1)

def makemovie(nsteps=1000,isteps=2):
  if top.it == 0: visualize(lwrite=1)
  while top.it < nsteps:
    step(isteps)
    visualize(lwrite=1)


makemovie()
exit(0)
