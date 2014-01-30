from warp import *
from em3dsolver import *
from tunnel_ionization import *
from boosted_frame import *
from species import *
import PRpickle as PR
import PWpickle as PW
import os
home=os.getenv('HOME')
# --- flags turning off unnecessary diagnostics (ignore for now)
top.ifzmmnt = 0
top.itmomnts = 0
top.itplps = 0
top.itplfreq = 0
top.zzmomnts = 0
top.zzplps = 0
top.zzplfreq = 0
top.nhist = top.nt
top.iflabwn = 0
w3d.lrhodia3d = false
w3d.lgetese3d = false
w3d.lgtlchg3d = false
#EnableAll()


#-------------------------------------------------------------------------------
# main parameters
#-------------------------------------------------------------------------------
#dim = "3d"                 # 3D calculation
dim = "2d"                 # 2D calculation 
#dim = "1d"                 # 1D calculation 
dpi=100                     # graphics resolution
l_test             = 0      # Will open output window on screen 
                            # and stop before entering main loop.
l_gist             = 1      # Turns gist plotting on/off
l_restart          = false  # To restart simulation from an old run (works?)
restart_dump       = ""     # dump file to restart from (works?)
l_moving_window    = 1      # on/off (Galilean) moving window
l_plasma           = 1      # on/off plasma
l_ions             = 1      # on/off plasma ions
l_external_field   = 1      # on/off external field for the injection pulse
l_beam             = 0      # on/off electron beam
l_injectplane      = 1      # on/off beam injection through plane
l_usesavedist      = 0      # if on, uses dump of beam particles distribution
savedist           = home+'/runs/warp/lhc/quasistatic/unitdist4sym300000' # beam initial distribution file
svstride           = 100    # loads only every svstride particles from dump
l_smooth           = 1      # on/off smoothing of current density
l_laser            = 1      # on/off laser
l_pdump            = 0      # on/off regular dump of beam data
stencil            = 0      # 0 = Yee; 1 = Yee-enlarged (Karkkainen) on EF,B; 2 = Yee-enlarged (Karkkainen) on E,F 
                            # use 0 or 1; 2 does not verify Gauss Law
if dim=="1d":stencil=0
dtcoef             = 0.25/4     # coefficient to multiply default time step that is set at the EM solver CFL
top.depos_order    = 3      # particles deposition order (1=linear, 2=quadratic, 3=cubic)
top.efetch         = 4      # field gather type (1=from nodes "momentum conserving"; 4=from Yee mesh "energy conserving")
nzstations         = 50     # number of beam diag z-stations
l_pselect          = 0      # on/off selection of particles (i.e. remove halo) for diagnostics
top.runid          = "lpa_basic"                         # run name
top.pline1         = "basic lpa"                         # comment line on plots
top.runmaker       = "J.-L. Vay,"                        # run makers
top.lrelativ       = true                                # on/off relativity (for particles push)
top.pgroup.lebcancel_pusher=true                         # flag for particle pusher (0=Boris pusher; 1=Vay PoP 08 pusher)
l_verbose          = 0                                   # verbosity level (0=off; 1=on)

#-------------------------------------------------------------------------------
# diagnostics parameters + a few other settings
#-------------------------------------------------------------------------------
hist_freq          = 50    # frequency (in time steps) of beam history data saving
live_plot_freq     = 2800*2  #2400  # frequency (in time steps) of live plots (off is l_test is off)

#-------------------------------------------------------------------------------
# boosted frame
#-------------------------------------------------------------------------------
gammafrm           = 1.
betafrm            = sqrt(1.-1./gammafrm**2)
if gammafrm>1.: # turns ON plasma ions if in boosted frame
  l_ions = 1
  l_moving_window = 1

#-------------------------------------------------------------------------------
# some units for convenience (& clarity)
#-------------------------------------------------------------------------------
microns            = 1.e-6
femtoseconds       = 1.e-15

lambda_laser_lab   = 5.0*microns             # wavelength 
lambda_laser_lab1  = 0.4*microns             # wavelength for the second laser



#-------------------------------------------------------------------------------
# plasma density, length and ramps 
#-------------------------------------------------------------------------------
# --- in lab frame
dfact             = 1.                                     # coefficient factor for plasma density (for scaled simulations)
dens0lab          = dfact*2.0e23                           # plasma density (flat section)
wplab             = sqrt(dens0lab*echarge**2/(eps0*emass)) # plasma frequency
kplab             = wplab/clight                           # plasma wavenumber
lambda_plasma_lab = 2.*pi/kplab                            # plasma wavelength
length_pramp_lab  = 15.*lambda_laser_lab/(dfact**1.5) #15.                     # plasma entrance ramp length
length_pramp_exit_lab = 15.*lambda_laser_lab/(dfact**1.5) #15.                  # plasma exit ramp length
Lplasma_lab       = 10000.*lambda_laser_lab/dfact**1.5 #500.                        # plasma total length (including ramps)
zstart_plasma_lab = 30.*lambda_laser_lab/dfact**1.5 #40.
# --- in boosted frame
dens0             = dens0lab*gammafrm                      # plasma density
length_pramp      = length_pramp_lab/gammafrm              # plasma ramp length
length_pramp_exit = length_pramp_exit_lab/gammafrm         # plasma ramp length
Lplasma           = Lplasma_lab/gammafrm                   # plasma total length
zstart_plasma     = zstart_plasma_lab/gammafrm             # plasma starting point

#-------------------------------------------------------------------------------
# ions density, length and ramps 
#-------------------------------------------------------------------------------
# --- in lab frame
dens_ions_lab     = 0.1*dens0lab                          # ions density (flat section)
length_iramp_lab  = 5.*lambda_laser_lab/(dfact**1.5) #4.                      # ions entrance ramp length
length_iramp_exit_lab = 5.*lambda_laser_lab/(dfact**1.5)#4.                  # ions exit ramp length
Lions_lab         = 20.*lambda_laser_lab/dfact**1.5#16.                        # ions total length (including ramps)
zstart_ions_lab   = 60.*lambda_laser_lab/(dfact**1.5)#80.    
# --- in boosted frame
dens_ions         = dens_ions_lab*gammafrm                 # ions density
length_iramp      = length_iramp_lab/gammafrm              # ions ramp length
length_iramp_exit = length_iramp_exit_lab/gammafrm         # ions ramp length
Lions             = Lions_lab/gammafrm                     # ions total length
zstart_ions       = zstart_ions_lab/gammafrm               # ions starting point

#-------------------------------------------------------------------------------
# laser parameters
#-------------------------------------------------------------------------------
# --- in lab frame
KP_L               = 2.0                    # normalized length - standard gaussian form P=P0*exp(-2xi^2/L^2)  
KP_SIGMA           = 3.02                   # normalized transverse spot size - standard gaussian form I=I0*exp(-2r^2/SIGMA^2)  
laser_radius       = KP_SIGMA / kplab         
laser_waist        = laser_radius*2.354820   # radius -> FWHM
laser_length_lab   = 4.75*lambda_laser_lab #4.75    # laser length 
laser_duration_lab = laser_length_lab/clight # laser duration
laser_risetime_lab = 3.*laser_duration_lab   # 
laser_dura_fwhm_lab= 1.1774*laser_duration_lab
laser_polangle     = 0                    # polarization (0=aligned with x; pi/2=aligned with y)
a0                 = 1.17                     # normalized potential vector (amplitude)
k0lab              = 2.*pi/lambda_laser_lab
w0lab              = k0lab*clight
ZR                 = 0.5*k0lab*(laser_waist**2)   # Rayleigh length
zstart_laser_lab   = 0.

# --- in boosted frame
lambda_laser       = lambda_laser_lab*gammafrm*(1.+betafrm)   # wavelength 
laser_duration     = laser_duration_lab*gammafrm*(1.+betafrm)
laser_length       = laser_duration/clight
laser_risetime     = 3.*laser_duration
laser_dura_fwhm    = 1.1774*laser_duration
zstart_laser       = 0. 
k0                 = 2.*pi/lambda_laser
w0                 = k0*clight
Eamp               = a0*w0*emass*clight/echarge
Bamp               = Eamp/clight 
if l_laser==0:
  Eamp*=1.e-15
  Bamp*=1.e-15


#-------------------------add second laser---------------------------------
KP_L1               = 2.0                     # normalized length - standard gaussian form P=P0*exp(-2xi^2/L^2)                        
KP_SIGMA1           = 2.*pi/15.                     # normalized transverse spot size - standard gaussian form I=I0*exp(-2r^2/SIGMA^2)       
laser_radius1       = KP_SIGMA1 / kplab 
laser_waist1        = laser_radius1*2.354820   # radius -> FWHM                                                                         
laser_length_lab1   = 10.*lambda_laser_lab1     # laser length                                                                          
laser_duration_lab1 = laser_length_lab1/clight # laser duration (FWHM)                                                                  
laser_risetime_lab1 = 3.*laser_duration_lab1
laser_dura_fwhm_lab1= 1.1774*laser_duration_lab1
laser_polangle1     = 0                    # polarization (0=aligned with x; pi/2=aligned with y)                                   
a1                 =0.135                    # normalized potential vector (amplitude)                                                
k1lab              = 2.*pi/lambda_laser_lab1 
w1lab              = k1lab*clight 
ZR1                 = 0.5*k1lab*(laser_waist1**2)   # Rayleigh length
#zstart_laser_lab1  = (7.6-28.)*lambda_laser_lab                                                                   
# --- in boosted frame                                                                                                                
lambda_laser1       = lambda_laser_lab1*gammafrm*(1.+betafrm)   # wavelength                                                            
laser_duration1     = laser_duration_lab1*gammafrm*(1.+betafrm) 
laser_length1       = laser_duration1/clight
laser_risetime1     = 3.*laser_duration1
laser_dura_fwhm1    = 1.1774*laser_duration1
#zstart_laser1  = zstart_laser_lab1/gammafrm 
k1                 = 2.*pi/lambda_laser1 
w1                 = k1*clight 
Eamp1               = a1*w1*emass*clight/echarge 
Bamp1               = Eamp1/clight 
if l_laser==0: 
  Eamp1*=1.e-15 
  Bamp1*=1.e-15 


#-------------------------------------------------------------------------------
# plasma cont'd
#-------------------------------------------------------------------------------
K                   = k0lab/kplab
Ld=LINEAR_DEPHASING = 0.5*lambda_plasma_lab**3/lambda_laser_lab**2             # linear dephasing length
BETAG_LINEAR_LAB    = sqrt(1-(1./K)*(1./K))                                    # linear beta of wake in lab
GAMMAG_LINEAR_LAB   = 1./sqrt(1-BETAG_LINEAR_LAB*BETAG_LINEAR_LAB)             # linear gamma of wake in lab
BETAG_LINEAR        = (BETAG_LINEAR_LAB-betafrm)/(1.-BETAG_LINEAR_LAB*betafrm) # linear beta of wake in simulation frame
GAMMAG_LINEAR       = 1./sqrt(1-BETAG_LINEAR*BETAG_LINEAR)                     # linear gamma of wake in simulation frame
kp                  = kplab*(gammafrm*(1.-BETAG_LINEAR_LAB*betafrm))
lambda_plasma       = 2.*pi/kp
densc               = emass*eps0*w0lab**2/echarge**2                           # critical density

#-------------------------------------------------------------------------------
# print some plasma parameters to the screen
#-------------------------------------------------------------------------------
print("the laser group velocity is: ")
print BETAG_LINEAR*clight
print("the laser spot size is: ")
print laser_waist
print("the Rayleigh length is: ")
print ZR
print("the laser wavelength is: ")
print lambda_laser
print("the second laser spot size is: ")
print laser_waist1
print("the Rayleigh length is: ")
print ZR1
print("the laser wavelength is: ")
print lambda_laser1
print("the plasma wavelength is: ")
print lambda_plasma
print("the plasma length is: ")
print Lplasma

#-------------------------------------------------------------------------------
# e-beam
#-------------------------------------------------------------------------------
# --- in lab frame
E_BEAM_GAMMA      = GAMMAG_LINEAR_LAB*1.5
E_BEAM_ENERGY_MEV = 0.511*(E_BEAM_GAMMA-1.)
E_BEAM_BETA       = sqrt(1.- 1./(E_BEAM_GAMMA*E_BEAM_GAMMA))
E_BEAM_U          = E_BEAM_GAMMA * E_BEAM_BETA
E_BEAM_RADIUS     = 0.825e-6/sqrt(dfact)
E_BEAM_LENGTH     = 0.85e-6/sqrt(dfact)

# --- transverse spread (RMS Gaussian)
GAMMAVXSIGMA = 0.
GAMMAVYSIGMA = 0.
# --- longitudinal spread (RMS Gaussian)
GAMMAVZSIGMA = 0.

E_BEAM_DENSITY_PEAK = 1.0e10
E_BEAM_PHASE        = 5.*pi/4.
E_BEAM_DISTANCE_BEHIND_LASER = E_BEAM_PHASE/kplab

#-------------------------------------------------------------------------------
# number of grid cells
#-------------------------------------------------------------------------------
# --- transverse
nx = 1200
# --- longitudinal
nzplambda =60 

#-------------------------------------------------------------------------------
# number of plasma macro-particles/cell
#-------------------------------------------------------------------------------
nppcellx = 1#2
nppcelly = 1
nppcellz = 2

if dim=="2d":
  nppcelly = 1
if dim=="1d":
  nppcellx = nppcelly = 1

#-------------------------------------------------------------------------------
# grid dimensions, nb cells and BC
#-------------------------------------------------------------------------------
w3d.xmmax = 3.0*laser_radius
w3d.xmmin = -w3d.xmmax
w3d.ymmax = w3d.xmmax
w3d.ymmin = -w3d.ymmax
w3d.nx = nx
w3d.ny = w3d.nx    
w3d.dy = (w3d.ymmax-w3d.ymmin)/w3d.ny
if dim in ["1d"]:
    w3d.nx = 2
    w3d.xmmin = -float(w3d.nx)/2
    w3d.xmmax = float(w3d.nx)/2
if dim in ["1d","2d"]:
    w3d.ny = 2
    w3d.ymmin = -float(w3d.ny)/2
    w3d.ymmax = float(w3d.ny)/2
w3d.zmmin = -3.5*lambda_plasma 
w3d.zmmax = 3.*lambda_laser#/nzplambda
w3d.nz = int((w3d.zmmax-w3d.zmmin)*nzplambda/lambda_laser) 
w3d.dx = (w3d.xmmax-w3d.xmmin)/w3d.nx
w3d.dy = (w3d.ymmax-w3d.ymmin)/w3d.ny
w3d.dz = (w3d.zmmax-w3d.zmmin)/w3d.nz

# --- enforces dx = dz (= dy) in boosted frame
if 0:#gammafrm>1:
  if w3d.dz<w3d.dx:
    # x (and y) dimensions are rescaled if needed
    w3d.dx = w3d.dz
    w3d.nx = 2*nint(0.5*(w3d.xmmax-w3d.xmmin)/w3d.dx)
    w3d.xmmax = w3d.nx*w3d.dx/2
    w3d.xmmin = -w3d.xmmax
    if dim=='3d':
      w3d.dy = w3d.dz
      w3d.ny = 2*nint(0.5*(w3d.ymmax-w3d.ymmin)/w3d.dy)
      w3d.ymmax = w3d.ny*w3d.dy/2
      w3d.ymmin = -w3d.ymmax
  elif w3d.dx<w3d.dz:
    # z dimensions are rescaled if needed
    w3d.dz = w3d.dx
    w3d.nz = nint((w3d.zmmax-w3d.zmmin)/w3d.dz)
    w3d.zmmax = w3d.zmmin+w3d.nz*w3d.dz

if gammafrm>10.:
  if dim=="1d":
    dt0 = w3d.dz/clight
  if dim=="2d":
    dt0 = 1./(clight*sqrt(1./w3d.dx**2+1./w3d.dz**2))
  if dim=="3d":
    dt0 = 1./(clight*sqrt(1./w3d.dx**2+1./w3d.dy**2+1./w3d.dz**2))
  mydt = w3d.dz/(sqrt(2.)*clight)
  dtcoef = min(1.,mydt/dt0)

# --- sets field boundary conditions
w3d.bound0  = w3d.boundnz = openbc
#w3d.bound0  = w3d.boundnz = -1 # reflective longitudinal BC
w3d.boundxy = periodic   #-1                # reflective transverse BC

# --- sets particles boundary conditions
top.pboundxy = periodic      #reflect#absorb
top.pbound0  = absorb
top.pboundnz = absorb

if dim=="1d":
  w3d.boundxy = periodic
  top.pboundxy = periodic

#-------------------------------------------------------------------------------
# Plasma channel
#-------------------------------------------------------------------------------
#  Matched spot size in channel (a0)
#  Note - to decrease channel dN at SIGMA, make matchspot > SIGMA
#  then dN/dN_matched = (SIGMA^4/matchspot^4), e.g. 5.7% increase in matchspot -> 20% detuning of dN 
#  tor: matchspot/sigma = (dN_matched/dN)^0.25

channel_region=3.0*laser_radius
CHANNELFAC = 0.4 # factor by which density rise is less than matched density (to compensate self guide)
matchspot  = laser_radius*(1/CHANNELFAC)**0.25 # channel
max_parab_radius   = 0.8*channel_region
max_radius        = 0.96*channel_region
diff_density      = gammafrm*1.13e14/(matchspot**4)  #formula for the channel density without the radius squared term
if dim=="1d":diff_density=0.
norm_diff_density  = diff_density/dens0        #the normalized channel density without the radius squared term
max_diff_density  = norm_diff_density*(max_parab_radius**2) #the maximum normalized channel density

#-------------------------------------------------------------------------------
# set max time
#-------------------------------------------------------------------------------
tmaxlab = (Lplasma_lab)/clight
tmax = tmaxlab / gammafrm

#-------------------------------------------------------------------------------
# dump/plots intervals (in pico-seconds)
#-------------------------------------------------------------------------------
dump_intervals = tmax/20
beamdump_intervals = tmax/10
plot_intervals = tmax/10

#-------------------------------------------------------------------------------
# set graphics
#-------------------------------------------------------------------------------
if l_gist:
 if l_test:
  window(0,dpi=dpi)
 else:
  setup()
else:
 setup()
 
#-------------------------------------------------------------------------------
# set particles 
#-------------------------------------------------------------------------------
weight     = 1.*dens0*w3d.dx*w3d.dy*w3d.dz/(nppcellx*nppcelly*nppcellz) # weight of plasma macro-particles
weightbeam = 0. # needs to be fixed

# --- create e- beam species
if l_beam:
  beam = Species(type=Electron,weight=weightbeam)
# --- create plasma electron species
elec = Species(type=Electron,weight=weight)
# --- create plasma electron species
prot = Species(type=Proton,weight=weight)
# --- create plasma ion species
# --- in this example, initial charge state is +5 
# --- charge states +6 and +7 are also considered
ielec = Species(type=Electron,weight=weight)
if l_ions:
  nions = 2
  ions=[]
  for i in range(nions):
    ions.append(Species(type=Krypton,charge_state=i+8,weight=weight,name='ions %g'%i))

top.pgroup.ndts=nint(1./dtcoef)
ielec.ndts=1

#top.pgroup.sm*=1000000

#  ions.sm*=1.e20  # to give virtually infinite mass to ions
top.wpid = nextpid() # creates data space for variable weights (needed for plasma ramps)
top.depos_order[...] = top.depos_order[0,0] # sets deposition order of all species = those of species 0
top.efetch[...] = top.efetch[0] # same for field gathering
if dim in ["1d","2d"]:
  top.depos_order[1,:]=1
if dim=="1d":
  top.depos_order[0,:]=1

#-------------------------------------------------------------------------------
# set smoothing of current density
#-------------------------------------------------------------------------------
if l_smooth:
  # --- 1 time nilinear (0.25,0.5,0.25) + 1 time relocalization (-1, 3/2,-1.)
  npass_smooth = [[ 1 , 1 ],[ 0 , 0 ],[ 4 , 1 ]]
  alpha_smooth = [[ 0.5, 3./2],[ 0.5, 3.],[0.5, 3./1]]
  stride_smooth = [[ 1 , 1 ],[ 1 , 1 ],[ 1 , 1 ]]
  if dim=='1d':
    for i in range(len(npass_smooth[0])):
      npass_smooth[0][i]=0
  if dim in ['1d','2d']:
    for i in range(len(npass_smooth[0])):
      npass_smooth[1][i]=0
else:
  npass_smooth = [[ 0 ],[ 0 ],[ 0 ]]
  alpha_smooth = [[ 1.],[ 1.],[ 1.]]
  stride_smooth = [[ 1 ],[ 1 ],[ 1 ]]

#-------------------------------------------------------------------------------
# initializes WARP
#-------------------------------------------------------------------------------
top.fstype = -1 # sets field solver to None (desactivates electrostatic solver)
package('w3d');generate()
#-------------------------------------------------------------------------------
# set a few shortcuts
#-------------------------------------------------------------------------------
pg = top.pgroup

#-------------------------------------------------------------------------------
# sets tunnel ionizations
#-------------------------------------------------------------------------------
tunnel_ioniz = TunnelIonization(stride=1)
for i in range(len(ions)-1):
  # --- For each ionization level, add tunnel ionization event by registering 
  # --- incident and emitted species.
  tunnel_ioniz.add(incident_species = ions[i],
                   emitted_species  = [ ions[i+1] , ielec ])

#-------------------------------------------------------------------------------
# set input plasma density arrays
#-------------------------------------------------------------------------------
# --- set intial positions and weights
zpos = zstart0 = 0.
nppcell = nppcellx*nppcelly*nppcellz
dx = w3d.dx/nppcellx
dy = w3d.dy/nppcelly
dz = w3d.dz/nppcellz
nx = nppcellx*w3d.nx
ny = nppcelly*w3d.ny
nz = nppcellz

if dim in ["1d","2d"]:
  if dim=="1d":
    nx=1
    xp0,zp0 = getmesh2d(0.,dx,0,
                        dz/2,dz,nz-1)
    yp0 = xp0*0.
  else:
    xp0,zp0 = getmesh2d(w3d.xmmin+dx/2,dx,nx-1,
                        -w3d.dz+dz/2,dz,nz-1)
    yp0 = xp0*0.
else:
  xp0,yp0,zp0 = getmesh3d(w3d.xmmin+dx/2,dx,nx-1,
                       w3d.ymmin+dy/2,dy,ny-1,
                       -w3d.dz+dz/2,dz,nz-1)

zp0-=minnd(zp0) # ensures that zp0 starts at 0

# --- transform to 1D arrays
xp0=xp0.flatten()
yp0=yp0.flatten()
zp0=zp0.flatten()

# --- select particles within computational box of local processor
ii=compress((xp0>=w3d.xmminlocal) & (xp0<w3d.xmmaxlocal) & \
            (yp0>=w3d.ymminlocal) & (yp0<w3d.ymmaxlocal),arange(len(xp0)))
xp0=take(xp0,ii)
yp0=take(yp0,ii)
zp0=take(zp0,ii)

# --- select particles within max_radius
rp0=sqrt(xp0**2+yp0**2)
ii=compress(rp0<=max_radius,arange(len(rp0)))
xp0=take(xp0,ii)
yp0=take(yp0,ii)
zp0=take(zp0,ii)
rp0=sqrt(xp0**2+yp0**2)

# --- set the transverse profile
def plasma_trans_profile(r):
    wp = ones(shape(r)[0])
    slope = -(1+max_diff_density)/(max_radius-max_parab_radius)
    intercept = -max_radius*slope
    wp = where(r<max_parab_radius,1+norm_diff_density*(r**2),slope*r+intercept)
    wp = where(wp<0.,0.,wp)
    return wp

wp0=plasma_trans_profile(rp0)

def plasma_long_profile(wp0,z,zstart):
  zp=z-zstart
  wp = where((zp<length_pramp) & (zp>0.),wp0*zp/length_pramp,wp0) # plasma density rises as linear ramp
  wp = where((zp>Lplasma-length_pramp_exit) & (zp<Lplasma),wp0*(1.-(zp-Lplasma+length_pramp_exit)/length_pramp_exit),wp) # plasma density falls as linear ramp
  wp = where((zp<=0.) | (zp>=Lplasma),0.,wp)
  return wp

def ions_long_profile(wp0,z,zstart):
  zp=z-zstart
  wp = where((zp<length_iramp) & (zp>0),wp0*zp/length_iramp,wp0) # plasma density rises as linear ramp
  #wp = where((zp>=length_iramp) & (zp<=Lions-length_iramp_exit),wp0,wp0)
  wp = where((zp>Lions-length_iramp_exit) & (zp<Lions),wp0*(1.-(zp-Lions+length_iramp_exit)/length_iramp_exit),wp) # plasma density falls as linear ramp
  wp = where((zp<=0.) | (zp>=Lions),0.,wp)
  return wp*dens_ions/dens0

def pldens():
  # --- plot longitudinal and transverse plasma profile
  nz = 1000
  za = arange(-0.,Lplasma,(Lplasma)/nz)
  we = ones(shape(za)[0])
  wp = ones(shape(za)[0])
  wi = ones(shape(za)[0])
  for i,z in enumerate(za.tolist()):
    we[i] = plasma_long_profile(1.,z,zstart_plasma)
    wi[i] = ions_long_profile(1.,z,zstart_ions)
    wp[i] = we[i]-wi[i]*(ions[0].charge/echarge)
  plsys(9)
  pla(wp,za*gammafrm*1000,width=3,color=blue)
  pla(wi,za*gammafrm*1000,width=3,color=green)
  pla(we,za*gammafrm*1000,width=3,color=red,type='dash')
  limits(0,za[-1]*gammafrm*1000,0,1.1)
  ptitles('longitudinal density profile','z (mm)','')
  plsys(10)
  r = arange(0,w3d.xmmax,w3d.xmmax/nz)
  wp = plasma_trans_profile(r)
  pla(wp,r*1.e6,width=3,color=red)
  limits(0,max_radius*1.e6,0,1.1*max(wp))
  ptitles('radial density profile','r (microns)','')

# --- defines subroutine injecting plasma
def loadplasma():
 global zstart0,zpos,xp0,yp0,zp0,wp0,l_moving_window
 while(zpos>=zstart0 and (zpos+betafrm*clight*top.time<Lplasma)): 
    z0 = zstart0+zp0
    # --- sets ramp by adjusting weight
    zi = z0+betafrm*clight*top.time
    # --- get electrons weight factor
    we = plasma_long_profile(wp0,zi,zstart_plasma)
    # --- get ions weight factor
    wi = ions_long_profile(wp0,zi,zstart_ions)
    # --- set protons weight factor
    wp = we-wi*(ions[0].charge/echarge)
    if dim<>'1d':wp*=plasma_trans_profile(rp0)
    # --- sets velocity
    vx = 0.#001*clight*ranf(dz)
    vy = 0.#001*clight*ranf(dz)
    vz = -betafrm*clight#+0.001*clight*ranf(dz)
    # --- sets positions
    x = xp0.copy()
    y = yp0.copy()
    z = z0.copy()
    if any(wp)>0.:
      # --- inject electrons
      elec.addpart(x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,w=we,lallindomain=false)
      prot.addpart(x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,w=wp,lallindomain=false)
      if any(wi)>0.:
        # --- inject ions at same locations
        ions[0].addpart(x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,w=wi,lallindomain=false)
    ladd=1
    zstart0+=w3d.dz
 if l_moving_window:
#    zpos+=top.vbeamfrm*top.dt
    zpos+=clight*top.dt # was not tested for a long time; is this correct?
 else:
    zpos+=clight*top.dt # was not tested for a long time; is this correct?
 zstart0-=betafrm*clight*top.dt

if l_plasma:installuserinjection(loadplasma)



#-------------------------------------------------------------------------------
# set laser pulse shape
#-------------------------------------------------------------------------------
def laser_amplitude(time):
 global laser_total_length,Eamp,laser_duration,laser_risetime
 l = laser_duration
 l_rise=laser_risetime
 return Eamp*exp(-((time-l_rise)/l)**2)

def laser_profile(x,y):
  global laser_waist
  r2 = x**2 + y**2
  rw = laser_radius  
  return exp(-r2/rw**2)

#-------------------------------------------------------------------------------
# set laser amplitude by combining the pulse shape, laser profile, and laser phase
#-------------------------------------------------------------------------------
def laser_func(x,y,t):
  # --- returns components [Ex,Ey] versus x,y,t.
  global laser_amplitude,laser_phase,laser_profile,k0,w0
  em.laser_source_z=zstart_laser
  E = laser_amplitude(t)*laser_profile(x,y)
  angle = w0*t   #k0*x-w0*t   #2.*pi*t*clight/lambda_laser
  return [E*cos(angle),E*sin(angle)]   # for circularly polarized laser
  #return [0,E*sin(angle)]               # for linearly polarized laser
  


#------------------------------ add second laser--------------------------------#----------------set the time delay between two lasers--------------------------
T_generate_laser=laser_risetime-laser_risetime1#2*laser_dura_fwhm  # time for generate laser pulses
Tstart_laser1=T_generate_laser
zstart_laser1=zstart_laser-21.25*lambda_laser
vg=sqrt(1-kplab*kplab/k1lab/k1lab)*clight


def laser_amplitude1(time):
 global laser_total_length1,Eamp1,laser_duration1,laser_risetime1
 l1=laser_duration1
 l_rise1=laser_risetime1
 return Eamp1*exp(-((time-l_rise1)/l1)**2)

def laser_profile1(x,y):
 global laser_waist1
 r21=x**2+y**2
 rw1=laser_radius1
 return exp(-r21/rw1**2)

def laser_func1(x,y,t):
 global laser_amplitude1,laser_phase1,laser_profile1,k1,w1,ZR1,vg,zstart_ions_lab
 em.laser_source_z=zstart_laser1
 z=vg*(top.time-laser_risetime)+zstart_laser1
 if z<=zstart_ions_lab:z=0
 else:z=z-zstart_ions_lab
 #print vg
 #print clight
 diffrac=1./sqrt(1+z*z/ZR1/ZR1)
 diffrac1=diffrac*diffrac
# E1=laser_amplitude1(t-Tstart_laser1)*laser_profile1(x,y)
 E1=laser_amplitude1(t-Tstart_laser1)*diffrac*pow(laser_profile1(x,y),diffrac1)
 angle1=w1*(t-Tstart_laser1)  #k1*zstart_laser1-w1*(t-Tstart_laser1) #2.*pi*(t-Tstart_laser1)*clight/lambda_laser1
 #return [E1*cos(angle1),E1*sin(angle1)]
 return [0,E1*sin(angle1)] #Ey
 #return [E1*sin(angle1),0] #Ex
 

#-------------------------------------------------------------------------------
# initializes main field solver block
#-------------------------------------------------------------------------------
laser_func_dict={}
laser_func_dict[0]=laser_func
if not l_external_field:
  laser_func_dict[1]=laser_func1

em = EM3D(       laser_func=laser_func_dict,    # laser_func=laser_func,
                 laser_source_z=0.,
                 laser_source_v=-betafrm*clight,
                 laser_polangle=laser_polangle,
                 laser_mode=2,
                 laser_emax=10*Eamp,
                 stencil=stencil,
                 npass_smooth=npass_smooth,
                 alpha_smooth=alpha_smooth,
                 stride_smooth=stride_smooth,
                 l_2dxz=dim=="2d",
                 l_1dz=dim=="1d",
                 dtcoef=dtcoef,
#                 l_enableovercycle=True,
                 l_verbose=l_verbose)

def add_external_laser():
    global vg
    #vg=sqrt(1-kplab*kplab/k1lab/k1lab)*clight
    il = w3d.jmin
    iu = w3d.jmax
    if iu<=il:return
    x = top.pgroup.xp[il:iu]
   # print me,il,iu,top.pgroup.yp[il:iu]
    y = top.pgroup.yp[il:iu]
    z = top.pgroup.zp[il:iu]
    t = top.time-(z-zstart_laser1)/vg
    ex1,ey1=laser_func1(x,y,t)
    bx1=-ey1/clight
    by1=ex1/clight
    top.pgroup.ex[il:iu] += ex1
    top.pgroup.ey[il:iu] += ey1
    top.pgroup.bx[il:iu] += bx1 
    top.pgroup.by[il:iu] += by1

if l_external_field:
  installothereuser(add_external_laser)

#-------------------------------------------------------------------------------
# restarts from dump file
#-------------------------------------------------------------------------------
if l_restart:
  restore(dump_file)

# --- load diagnostics
execfile('lpa_basic_diags.py')

#-------------------------------------------------------------------------------
# intializes e- beam
#-------------------------------------------------------------------------------
if l_beam:
  # --- add beam particles
  np_beam = 4000
  top.vbeam = E_BEAM_BETA*clight
  if me==0: # --- do only if processor 0
   if np_beam==1:
    # --- loads single test electron
    beam.addpart(x=array([0.]),
                 y=array([0.]),
                 z=array([0.]),
                 vx=array([0.]),
                 vy=array([0.]),
                 vz=array([E_BEAM_GAMMA*E_BEAM_BETA*clight]),
                 gi=array([1./E_BEAM_GAMMA]),
                         lmomentum=True,
                         lallindomain=True)
   else:
    # --- loads e- beam
     if l_usesavedist:
       # --- loads distribution from file
       try:
         ff=PR.PR(savedist+'.pdb')
       except:
         ff=PR.PR(savedist+'.pyp')
       ux = ff.xp[::svstride]*GAMMAVXSIGMA
       uy = ff.yp[::svstride]*GAMMAVYSIGMA
       uz = ff.dp[::svstride]*GAMMAVZSIGMA+E_BEAM_GAMMA*E_BEAM_BETA*clight
       gi = 1./sqrt(1.+(ux**2+uy**2+uz**2)/clight**2)
       beam.addpart(ff.x[::svstride]*E_BEAM_RADIUS*2,
                    ff.y[::svstride]*E_BEAM_RADIUS*2,
                    ff.z[::svstride]*E_BEAM_LENGTH,
                    ux,uy,uz,gi=gi,
                    lmomentum=True,
                    lallindomain=True)
     else:
       # --- loads gaussian electron beam
       beam.add_gaussian_dist(np=np_beam,
                              deltax=E_BEAM_RADIUS*2*1,
                              deltay=E_BEAM_RADIUS*2*1,
                              deltaz=E_BEAM_LENGTH,
                              vthx=GAMMAVXSIGMA*1,
                              vthy=GAMMAVYSIGMA*1,
                              vthz=GAMMAVZSIGMA,
                              zmean=0.,
                              vzmean=E_BEAM_GAMMA*E_BEAM_BETA*clight,
                              lmomentum=True,
                              lallindomain=True,
                              zdist='regular')
  np_beam = beam.getn()
  # --- sets e- beam macro-particles weight
  if dim=="1d":
    beam.sw = (E_BEAM_DENSITY_PEAK*((E_BEAM_RADIUS*2)*w3d.dx*w3d.dy*E_BEAM_LENGTH*(2.*pi)))/np_beam
  if dim=="2d":
    beam.sw = (E_BEAM_DENSITY_PEAK*((E_BEAM_RADIUS*2)*w3d.dy*E_BEAM_LENGTH*(2.*pi)))/np_beam
  if dim=="3d":
    beam.sw = (E_BEAM_DENSITY_PEAK*((E_BEAM_RADIUS*2)**2*E_BEAM_LENGTH*(2.*pi)*1.5))/np_beam
  # --- install e- beam diagnostic routine
  if sum(pg.nps)>0:
    # --- set beam position in lab frame
    pg.zp += zstart0-E_BEAM_DISTANCE_BEHIND_LASER-0.5*laser_total_length_lab
  # --- transform particle positions and velocity to boosted frame
  bf=Boosted_Frame(gammafrm,l_setselfb=0)
  zinit = getz().copy()
  bf.boost(beam,l_inject_plane=l_injectplane,lallindomain=1,l_rmzmean=0,zinject=-5.*w3d.dz)
  particleboundaries3d(top.pgroup,-1,False)

#-------------------------------------------------------------------------------
# register solver
#-------------------------------------------------------------------------------
print 'register solver'
registersolver(em)
print 'done'

#-------------------------------------------------------------------------------
# sets moving window velocity
#-------------------------------------------------------------------------------
if l_moving_window:
  #top.vbeamfrm=BETAG_LINEAR*clight
  top.vbeamfrm=clight
#-------------------------------------------------------------------------------
# set a few shortcuts
#-------------------------------------------------------------------------------
el=elec             
f = em.fields
#bhist=ones([nions,w3d.nz+1])
#bhist[1:,:]=0.
#def accuhist():
 # global bhist
 # for i in range(nions-2,-1,-1):
  #  Ex = em.gatherex()
  #  Ey = em.gatherey()
   # if me==0:
    #  E = sqrt(Ex*Ex+Ey*Ey)
     # E0 = where(E<=1.e-10*Eamp,1.e-10*Eamp,E)
     # dn = where(E==0.,0.,top.dt*tunnel_ioniz.GetADKrateSI(abs(E0),i+8,ions[i].type,top.dt))
     # bhist[i+1,:]+=dn*bhist[i]
     # bhist[i,:]-=dn*bhist[i]

#installafterstep(accuhist)


window(1,hcp='lineout.cgm',dump=1,display='')
window(2,hcp='Ex.cgm',dump=1,display='') 


def liveplots():
  if top.it%live_plot_freq==0:
     #save_field_data()
     save_emittance()
     window(0)
     if l_test:fma()
    #nz=500
    #d = (ions[1].get_density(nz=nz)+ions[2].get_density(nz=nz))/dens0
    #d=ions[1].get_density(nz=nz)/dens_ions
    #dz=(w3d.zmmax-w3d.zmmin)/nz
    #z=w3d.zmmin+arange(nz+1)*dz
    #plsys(3);pla(d,z/lambda_laser,color=blue)
    #z=w3d.zmmin+arange(w3d.nz+1)*w3d.dz
    #pla((bhist[1,:]),z/lambda_laser,color=red)
    #pla(d,z/lambda_laser,color=blue)
    #pzxey(view=3,msize=1)
    #pzxex(view=4,msize=1)
    #pzxez(view=5,msize=1)
    #plke(view=6)
    #plez(10)
     if dim=="1d":
       plex(14)
       pldens1d(15)
       plzuz(16)
     else:
 #      plex(14)
 #      pldenlineout(15)
       plex2d(3)
#       plez(4)
       plex2dwake(5)
  #     pldens(6)
       plzuz(6)
     if not l_test:fma()
     window(1)
#     plex(3)
     pldenlineout(4)
     if dim<>"1d":
       pldens(5)
     plspectrum(6)
     if not l_test:fma()
     window(2)
     plex_x_z(3)
     plx_z_ex(4)
     plxux(5)
    # plx_z_ex_p(6)
     if l_test:
       refresh()
     else:
       fma()
   
installafterstep(liveplots)

print '\Initialization complete\n'
step(3001*280)

#step(nint(560/dtcoef*float(nzplambda)/16))

