import os
import sys
import argparse
import math
import numpy
import warpoptions
warpoptions.lskipoptions = 1
import warp
import em3dsolver
import boosted_frame
import species
import PRpickle as PR
import PWpickle as PW
from WarpVisItSimulation import WarpVisItSimulation
from WarpVisItFilteredSpecies import WarpVisItSpeciesFilters

#----------------------------------------------------------------------------
def NewWarpVisItSimulation(args=[]):
    """
    factory that creates an obect implementing WarpVisItSimulation
    class interface
    """
    return LPATwoColorSimulation(args)

#----------------------------------------------------------------------------
class LPATwoColorSimulation(WarpVisItSimulation):
    """
    Implementation of the WarpVisItSimulation class for two color
    lpa sim example.
    """
    def __init__(self,args=[]):
        global bigRun
        # parse command line args
        self.bigRun = False

        # parse command line args
        ap = argparse.ArgumentParser(usage=argparse.SUPPRESS,prog='LPATwoColorSimulation',add_help=False)
        ap.add_argument('--big-run',default=self.bigRun,action='store_true')
        opts = vars(ap.parse_known_args(args)[0])
        self.bigRun = opts['big_run']
        bigRun = self.bigRun

        # call base class constructor
        WarpVisItSimulation.__init__(self,args)

        # setup vis scripts
        self.AddRenderScript('4Views','render-four-views.py')
        return

    #------------------------------------------------------------------------
    def Initialize(self):
        """
        Setup IC and start the simulation.
        """
        # call scientist provided script to intialize warp
        initlpa()

        # call default implementation
        WarpVisItSimulation.Initialize(self)

        print 'initialization complete'
        return


# ---------------------------------------------------------------------------
# helper functions
# ---------------------------------------------------------------------------
#Declare global variables used in various setup functions. The variables are
#initalized in the Initalize function

sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import tunnel_ionization

zstart0=None
zpos=None
xp0=None
yp0=None
zp0=None
wp0=None
l_moving_window=None
laser_total_length=None
Eamp=None
laser_duration=None
laser_risetime=None
laser_waist=None
laser_amplitude=None
laser_phase=None
laser_profile=None
k0=None
w0=None
vg=None
max_diff_density=None
max_radius=None
max_parab_radius=None
norm_diff_density=None
betafrm=None
Lplasma=None
zstart_plasma=None
zstart_ions=None
length_pramp = None
length_pramp_exit=None
length_iramp=None
Lions=None
length_iramp_exit=None
dens_ions=None
dens0=None
ions=None
dim=None
rp0=None
zstart_laser=None
em=None
laser_radius=None
laser_radius1 = None
live_plot_freq = None
elec=None
prot=None
ions=None
zstart_laser1=None
zstart_ions_lab=None
Tstart_laser1=None
laser_amplitude1=None
laser_phase1=None
laser_profile1=None
k1=None
w1=None
ZR1=None
laser_total_length1=None
Eamp1=None
laser_duration1=None
laser_risetime1=None
bigRun = None

def plasma_trans_profile(r):
    global max_diff_density, max_radius, max_parab_radius, norm_diff_density
    wp = numpy.ones(numpy.shape(r)[0])
    slope = -(1+max_diff_density)/(max_radius-max_parab_radius)
    intercept = -max_radius*slope
    wp = numpy.where(r<max_parab_radius,1+norm_diff_density*(r**2),slope*r+intercept)
    wp = numpy.where(wp<0.,0.,wp)
    return wp

def plasma_long_profile(wp0,z,zstart):
    global length_pramp, Lplasma, length_pramp_exit
    zp=z-zstart
    wp = numpy.where((zp<length_pramp) & (zp>0.),wp0*zp/length_pramp,wp0) # plasma density rises as linear ramp
    wp = numpy.where((zp>Lplasma-length_pramp_exit) & (zp<Lplasma),wp0*(1.-(zp-Lplasma+length_pramp_exit)/length_pramp_exit),wp) # plasma density falls as linear ramp
    wp = numpy.where((zp<=0.) | (zp>=Lplasma),0.,wp)
    return wp

def ions_long_profile(wp0,z,zstart):
    global length_iramp, Lions, length_iramp_exit, dens_ions, dens0
    zp=z-zstart
    wp = numpy.where((zp<length_iramp) & (zp>0),wp0*zp/length_iramp,wp0) # plasma density rises as linear ramp
    #wp = numpy.where((zp>=length_iramp) & (zp<=Lions-length_iramp_exit),wp0,wp0)
    wp = numpy.where((zp>Lions-length_iramp_exit) & (zp<Lions),wp0*(1.-(zp-Lions+length_iramp_exit)/length_iramp_exit),wp) # plasma density falls as linear ramp
    wp = numpy.where((zp<=0.) | (zp>=Lions),0.,wp)
    return wp*dens_ions/dens0

def pldens():
  # --- plot longitudinal and transverse plasma profile
  nz = 1000
  za = numpy.arange(-0.,Lplasma,(Lplasma)/nz)
  we = numpy.ones(numpy.shape(za)[0])
  wp = numpy.ones(numpy.shape(za)[0])
  wi = numpy.ones(numpy.shape(za)[0])
  for i,z in enumerate(za.tolist()):
    we[i] = plasma_long_profile(1.,z,zstart_plasma)
    wi[i] = ions_long_profile(1.,z,zstart_ions)
    wp[i] = we[i]-wi[i]*(ions[0].charge/warp.echarge)
  plsys(9)
  pla(wp,za*gammafrm*1000,width=3,color=blue)
  pla(wi,za*gammafrm*1000,width=3,color=green)
  pla(we,za*gammafrm*1000,width=3,color=red,type='dash')
  limits(0,za[-1]*gammafrm*1000,0,1.1)
  ptitles('longitudinal density profile','z (mm)','')
  plsys(10)
  r = numpy.arange(0,warp.w3d.xmmax,warp.w3d.xmmax/nz)
  wp = plasma_trans_profile(r)
  pla(wp,r*1.e6,width=3,color=red)
  limits(0,max_radius*1.e6,0,1.1*max(wp))
  ptitles('radial density profile','r (microns)','')

# --- defines subroutine injecting plasma
def loadplasma():
  global zstart0,zpos,xp0,yp0,zp0,wp0,l_moving_window, betafrm, Lplasma, zstart_plasma, zstart_ions, ions, dim, rp0, elec , prot, ions
  while(zpos>=zstart0 and (zpos+betafrm*warp.clight*warp.top.time<Lplasma)):
     z0 = zstart0+zp0
     # --- sets ramp by adjusting weight
     zi = z0+betafrm*warp.clight*warp.top.time
     # --- get electrons weight factor
     we = plasma_long_profile(wp0,zi,zstart_plasma)
     # --- get ions weight factor
     wi = ions_long_profile(wp0,zi,zstart_ions)
     # --- set protons weight factor
     wp = we-wi*(ions[0].charge/warp.echarge)
     if dim<>'1d':wp*=plasma_trans_profile(rp0)
     # --- sets velocity
     vx = 0.#001*warp.clight*ranf(dz)
     vy = 0.#001*warp.clight*ranf(dz)
     vz = -betafrm*warp.clight#+0.001*warp.clight*ranf(dz)
     # --- sets positions
     x = xp0.copy()
     y = yp0.copy()
     z = z0.copy()
     if any(wp)>0.:
       # --- inject electrons
       elec.addpart(x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,w=we,lallindomain=False)
       prot.addpart(x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,w=wp,lallindomain=False)
       if any(wi)>0.:
         # --- inject ions at same locations
         ions[0].addpart(x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,w=wi,lallindomain=False)
     ladd=1
     zstart0+=warp.w3d.dz
  if l_moving_window:
 #    zpos+=top.vbeamfrm*top.dt
     zpos+=warp.clight*warp.top.dt # was not tested for a long time; is this correct?
  else:
     zpos+=warp.clight*warp.top.dt # was not tested for a long time; is this correct?
  zstart0-=betafrm*warp.clight*warp.top.dt

def laser_amplitude(time):
    global laser_total_length,Eamp,laser_duration,laser_risetime
    l = laser_duration
    l_rise=laser_risetime
    return Eamp*numpy.exp(-((time-l_rise)/l)**2)

def laser_profile(x,y):
    global laser_waist, laser_radius
    r2 = x**2 + y**2
    rw = laser_radius
    return numpy.exp(-r2/rw**2)

def laser_func(x,y,t):
    # --- returns components [Ex,Ey] versus x,y,t.
    global laser_amplitude,laser_phase,laser_profile,k0,w0, zstart_laser, em
    em.laser_source_z=zstart_laser
    E = laser_amplitude(t)*laser_profile(x,y)
    angle = w0*t   #k0*x-w0*t   #2.*math.pi*t*warp.clight/lambda_laser
    return [E*math.cos(angle),E*math.sin(angle)]   # for circularly polarized laser
    #return [0,E*sin(angle)]               # for linearly polarized laser

def laser_amplitude1(time):
    global laser_total_length1,Eamp1,laser_duration1,laser_risetime1
    l1=laser_duration1
    l_rise1=laser_risetime1
    return Eamp1*numpy.exp(-((time-l_rise1)/l1)**2)

def laser_profile1(x,y):
    global laser_waist1, laser_radius1
    r21=x**2+y**2
    rw1=laser_radius1
    return numpy.exp(-r21/rw1**2)

def laser_func1(x,y,t):
    global laser_amplitude1,laser_phase1,laser_profile1,k1,w1,ZR1,vg,zstart_ions_lab, Tstart_laser1
    global bigRun
    em.laser_source_z=zstart_laser1
    z=vg*(warp.top.time-laser_risetime)+zstart_laser1
    if z<=zstart_ions_lab:z=0
    else:z=z-zstart_ions_lab
    #print vg
    #print warp.clight
    diffrac=1./math.sqrt(1+z*z/ZR1/ZR1)
    diffrac1=diffrac*diffrac
    # E1=laser_amplitude1(t-Tstart_laser1)*laser_profile1(x,y)
    E1=laser_amplitude1(t-Tstart_laser1)*diffrac*pow(laser_profile1(x,y),diffrac1)
    angle1=w1*(t-Tstart_laser1)  #k1*zstart_laser1-w1*(t-Tstart_laser1) #2.*math.pi*(t-Tstart_laser1)*warp.clight/lambda_laser1
    #return [E1*cos(angle1),E1*sin(angle1)]
    if bigRun:
        return [E1*numpy.sin(angle1),0] #Ex
    else:
        return [0,E1*numpy.sin(angle1)] #Ey

def add_external_laser():
    global vg, zstart_laser1
    #vg=math.sqrt(1-kplab*kplab/k1lab/k1lab)*warp.clight
    il = warp.w3d.jmin
    iu = warp.w3d.jmax
    if iu<=il:return
    x = warp.top.pgroup.xp[il:iu]
   # print me,il,iu,warp.top.pgroup.yp[il:iu]
    y = warp.top.pgroup.yp[il:iu]
    z = warp.top.pgroup.zp[il:iu]
    t = warp.top.time-(z-zstart_laser1)/vg
    ex1,ey1=laser_func1(x,y,t)
    bx1=-ey1/warp.clight
    by1=ex1/warp.clight
    warp.top.pgroup.ex[il:iu] += ex1
    warp.top.pgroup.ey[il:iu] += ey1
    warp.top.pgroup.bx[il:iu] += bx1
    warp.top.pgroup.by[il:iu] += by1


#----------------------------------------------------------------------------
# main initialization function
#----------------------------------------------------------------------------
def initlpa():
    global zstart0,zpos,xp0,yp0,zp0,wp0,l_moving_window, diagnosticsScript, laser_total_length,Eamp,laser_duration,laser_risetime, laser_waist, laser_amplitude,laser_phase,laser_profile,k0,w0, vg, max_diff_density, max_radius, max_parab_radius, norm_diff_density, betafrm, Lplasma, zstart_plasma, zstart_ions, length_pramp, length_pramp_exit, length_iramp, Lions, length_iramp_exit, dens_ions, dens0, ions, rp0, zstart_laser, em, laser_radius, laser_radius1, live_plot_freq, elec , prot, ions, zstart_laser1, zstart_ions_lab, Tstart_laser1, laser_amplitude1, laser_phase1, laser_profile1, k1, w1, ZR1, laser_total_length1, Eamp1,laser_duration1, laser_risetime1
    global bigRun

    # --- flags turning off unnecessary diagnostics (ignore for now)
    warp.top.ifzmmnt = 0
    warp.top.itmomnts = 0
    warp.top.itplps = 0
    warp.top.itplfreq = 0
    warp.top.zzmomnts = 0
    warp.top.zzplps = 0
    warp.top.zzplfreq = 0
    warp.top.nhist = warp.top.nt
    warp.top.iflabwn = 0
    warp.w3d.lrhodia3d = False
    warp.w3d.lgetese3d = False
    warp.w3d.lgtlchg3d = False
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
    l_gist             = 0      # Turns gist plotting on/off
    l_restart          = False  # To restart simulation from an old run (works?)
    restart_dump       = ""     # dump file to restart from (works?)
    l_moving_window    = 1      # on/off (Galilean) moving window
    l_plasma           = 1      # on/off plasma
    l_ions             = 1      # on/off plasma ions
    l_external_field   = 1      # on/off external field for the injection pulse
    l_beam             = 0      # on/off electron beam
    l_injectplane      = 1      # on/off beam injection through plane
    l_usesavedist      = 0      # if on, uses dump of beam particles distribution
    savedist           = 'unitdist4sym300000' # beam initial distribution file
    svstride           = 100    # loads only every svstride particles from dump
    l_smooth           = 1      # on/off smoothing of current density
    l_laser            = 1      # on/off laser
    l_pdump            = 0      # on/off regular dump of beam data
    stencil            = 0      # 0 = Yee; 1 = Yee-enlarged (Karkkainen) on EF,B; 2 = Yee-enlarged (Karkkainen) on E,F
                                # use 0 or 1; 2 does not verify Gauss Law
    if dim=="1d":stencil=0
                                # coefficient to multiply default time step that is set at the EM solver CFL
    if bigRun:
        dtcoef         = 0.25/4.0
    else:
        dtcoef         = 1.0

    warp.top.depos_order    = 3      # particles deposition order (1=linear, 2=quadratic, 3=cubic)
    warp.top.efetch         = 4      # field gather type (1=from nodes "momentum conserving"; 4=from Yee mesh "energy conserving")
    nzstations         = 50     # number of beam diag z-stations
    l_pselect          = 0      # on/off selection of particles (i.e. remove halo) for diagnostics
    warp.top.runid          = "lpa_basic"                         # run name
    warp.top.pline1         = "basic lpa"                         # comment line on plots
    warp.top.runmaker       = "J.-L. Vay,"                        # run makers
    warp.top.lrelativ       = True                                # on/off relativity (for particles push)
    warp.top.pgroup.lebcancel_pusher=True                         # flag for particle pusher (0=Boris pusher; 1=Vay PoP 08 pusher)
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
    betafrm            = math.sqrt(1.-1./gammafrm**2)
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
    wplab             = math.sqrt(dens0lab*warp.echarge**2/(warp.eps0*warp.emass)) # plasma frequency
    kplab             = wplab/warp.clight                           # plasma wavenumber
    lambda_plasma_lab = 2.*math.pi/kplab                            # plasma wavelength
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
    KP_L               = 2.0                    # normalized length - standard gaussian form P=P0*numpy.exp(-2xi^2/L^2)
    KP_SIGMA           = 3.02                   # normalized transverse spot size - standard gaussian form I=I0*numpy.exp(-2r^2/SIGMA^2)
    laser_radius       = KP_SIGMA / kplab
    laser_waist        = laser_radius*2.354820   # radius -> FWHM
    laser_length_lab   = 4.75*lambda_laser_lab #4.75    # laser length
    laser_duration_lab = laser_length_lab/warp.clight # laser duration
    laser_risetime_lab = 3.*laser_duration_lab   #
    laser_dura_fwhm_lab= 1.1774*laser_duration_lab
    laser_polangle     = 0                    # polarization (0=aligned with x; math.pi/2=aligned with y)
    a0                 = 1.17                     # normalized potential vector (amplitude)
    k0lab              = 2.*math.pi/lambda_laser_lab
    w0lab              = k0lab*warp.clight
    ZR                 = 0.5*k0lab*(laser_waist**2)   # Rayleigh length
    zstart_laser_lab   = 0.

    # --- in boosted frame
    lambda_laser       = lambda_laser_lab*gammafrm*(1.+betafrm)   # wavelength
    laser_duration     = laser_duration_lab*gammafrm*(1.+betafrm)
    laser_length       = laser_duration/warp.clight
    laser_risetime     = 3.*laser_duration
    laser_dura_fwhm    = 1.1774*laser_duration
    zstart_laser       = 0.
    k0                 = 2.*math.pi/lambda_laser
    w0                 = k0*warp.clight
    Eamp               = a0*w0*warp.emass*warp.clight/warp.echarge
    Bamp               = Eamp/warp.clight
    if l_laser==0:
      Eamp*=1.e-15
      Bamp*=1.e-15


    #-------------------------add second laser---------------------------------
    KP_L1               = 2.0                     # normalized length - standard gaussian form P=P0*numpy.exp(-2xi^2/L^2)
    KP_SIGMA1           = 2.*math.pi/15.                     # normalized transverse spot size - standard gaussian form I=I0*numpy.exp(-2r^2/SIGMA^2)
    laser_radius1       = KP_SIGMA1 / kplab
    laser_waist1        = laser_radius1*2.354820   # radius -> FWHM
    laser_length_lab1   = 10.*lambda_laser_lab1     # laser length
    laser_duration_lab1 = laser_length_lab1/warp.clight # laser duration (FWHM)
    laser_risetime_lab1 = 3.*laser_duration_lab1
    laser_dura_fwhm_lab1= 1.1774*laser_duration_lab1
    laser_polangle1     = 0                    # polarization (0=aligned with x; math.pi/2=aligned with y)
    a1                 =0.135                    # normalized potential vector (amplitude)
    k1lab              = 2.*math.pi/lambda_laser_lab1
    w1lab              = k1lab*warp.clight
    ZR1                 = 0.5*k1lab*(laser_waist1**2)   # Rayleigh length
    #zstart_laser_lab1  = (7.6-28.)*lambda_laser_lab
    # --- in boosted frame
    lambda_laser1       = lambda_laser_lab1*gammafrm*(1.+betafrm)   # wavelength
    laser_duration1     = laser_duration_lab1*gammafrm*(1.+betafrm)
    laser_length1       = laser_duration1/warp.clight
    laser_risetime1     = 3.*laser_duration1
    laser_dura_fwhm1    = 1.1774*laser_duration1
    #zstart_laser1  = zstart_laser_lab1/gammafrm
    k1                 = 2.*math.pi/lambda_laser1
    w1                 = k1*warp.clight
    Eamp1               = a1*w1*warp.emass*warp.clight/warp.echarge
    Bamp1               = Eamp1/warp.clight
    if l_laser==0:
      Eamp1*=1.e-15
      Bamp1*=1.e-15


    #-------------------------------------------------------------------------------
    # plasma cont'd
    #-------------------------------------------------------------------------------
    K                   = k0lab/kplab
    Ld=LINEAR_DEPHASING = 0.5*lambda_plasma_lab**3/lambda_laser_lab**2             # linear dephasing length
    BETAG_LINEAR_LAB    = math.sqrt(1-(1./K)*(1./K))                                    # linear beta of wake in lab
    GAMMAG_LINEAR_LAB   = 1./math.sqrt(1-BETAG_LINEAR_LAB*BETAG_LINEAR_LAB)             # linear gamma of wake in lab
    BETAG_LINEAR        = (BETAG_LINEAR_LAB-betafrm)/(1.-BETAG_LINEAR_LAB*betafrm) # linear beta of wake in simulation frame
    GAMMAG_LINEAR       = 1./math.sqrt(1-BETAG_LINEAR*BETAG_LINEAR)                     # linear gamma of wake in simulation frame
    kp                  = kplab*(gammafrm*(1.-BETAG_LINEAR_LAB*betafrm))
    lambda_plasma       = 2.*math.pi/kp
    densc               = warp.emass*warp.eps0*w0lab**2/warp.echarge**2                           # critical density

    #-------------------------------------------------------------------------------
    # print some plasma parameters to the screen
    #-------------------------------------------------------------------------------
    print("the laser group velocity is: ")
    print BETAG_LINEAR*warp.clight
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
    E_BEAM_BETA       = math.sqrt(1.- 1./(E_BEAM_GAMMA*E_BEAM_GAMMA))
    E_BEAM_U          = E_BEAM_GAMMA * E_BEAM_BETA
    E_BEAM_RADIUS     = 0.825e-6/math.sqrt(dfact)
    E_BEAM_LENGTH     = 0.85e-6/math.sqrt(dfact)

    # --- transverse spread (RMS Gaussian)
    GAMMAVXSIGMA = 0.
    GAMMAVYSIGMA = 0.
    # --- longitudinal spread (RMS Gaussian)
    GAMMAVZSIGMA = 0.

    E_BEAM_DENSITY_PEAK = 1.0e10
    E_BEAM_PHASE        = 5.*math.pi/4.
    E_BEAM_DISTANCE_BEHIND_LASER = E_BEAM_PHASE/kplab

    #-------------------------------------------------------------------------------
    # number of grid cells
    #-------------------------------------------------------------------------------
    # --- transverse
    # --- longitudinal
    if bigRun:
        nx = 1200
        nzplambda = 60
    else:
        nx = 1200/6
        nzplambda = 60/6


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
    warp.w3d.xmmax = 3.0*laser_radius
    warp.w3d.xmmin = -warp.w3d.xmmax
    warp.w3d.ymmax = warp.w3d.xmmax
    warp.w3d.ymmin = -warp.w3d.ymmax
    warp.w3d.nx = nx
    warp.w3d.ny = warp.w3d.nx
    warp.w3d.dy = (warp.w3d.ymmax-warp.w3d.ymmin)/warp.w3d.ny
    if dim in ["1d"]:
        warp.w3d.nx = 2
        warp.w3d.xmmin = -float(warp.w3d.nx)/2
        warp.w3d.xmmax = float(warp.w3d.nx)/2
    if dim in ["1d","2d"]:
        warp.w3d.ny = 2
        warp.w3d.ymmin = -float(warp.w3d.ny)/2
        warp.w3d.ymmax = float(warp.w3d.ny)/2
    warp.w3d.zmmin = -3.5*lambda_plasma
    warp.w3d.zmmax = 3.*lambda_laser#/nzplambda
    warp.w3d.nz = int((warp.w3d.zmmax-warp.w3d.zmmin)*nzplambda/lambda_laser)
    warp.w3d.dx = (warp.w3d.xmmax-warp.w3d.xmmin)/warp.w3d.nx
    warp.w3d.dy = (warp.w3d.ymmax-warp.w3d.ymmin)/warp.w3d.ny
    warp.w3d.dz = (warp.w3d.zmmax-warp.w3d.zmmin)/warp.w3d.nz

    # --- enforces dx = dz (= dy) in boosted frame
    if 0:#gammafrm>1:
      if warp.w3d.dz<warp.w3d.dx:
        # x (and y) dimensions are rescaled if needed
        warp.w3d.dx = warp.w3d.dz
        warp.w3d.nx = 2*warp.nint(0.5*(warp.w3d.xmmax-warp.w3d.xmmin)/warp.w3d.dx)
        warp.w3d.xmmax = warp.w3d.nx*warp.w3d.dx/2
        warp.w3d.xmmin = -warp.w3d.xmmax
        if dim=='3d':
          warp.w3d.dy = warp.w3d.dz
          warp.w3d.ny = 2*warp.nint(0.5*(warp.w3d.ymmax-warp.w3d.ymmin)/warp.w3d.dy)
          warp.w3d.ymmax = warp.w3d.ny*warp.w3d.dy/2
          warp.w3d.ymmin = -warp.w3d.ymmax
      elif warp.w3d.dx<warp.w3d.dz:
        # z dimensions are rescaled if needed
        warp.w3d.dz = warp.w3d.dx
        warp.w3d.nz = warp.nint((warp.w3d.zmmax-warp.w3d.zmmin)/warp.w3d.dz)
        warp.w3d.zmmax = warp.w3d.zmmin+warp.w3d.nz*warp.w3d.dz

    if gammafrm>10.:
      if dim=="1d":
        dt0 = warp.w3d.dz/warp.clight
      if dim=="2d":
        dt0 = 1./(warp.clight*math.sqrt(1./warp.w3d.dx**2+1./warp.w3d.dz**2))
      if dim=="3d":
        dt0 = 1./(warp.clight*math.sqrt(1./warp.w3d.dx**2+1./warp.w3d.dy**2+1./warp.w3d.dz**2))
      mydt = warp.w3d.dz/(math.sqrt(2.)*warp.clight)
      dtcoef = min(1.,mydt/dt0)

    # --- sets field boundary conditions
    warp.w3d.bound0  = warp.w3d.boundnz = warp.openbc
    #warp.w3d.bound0  = warp.w3d.boundnz = -1 # reflective longitudinal BC
    warp.w3d.boundxy = warp.periodic   #-1                # reflective transverse BC

    # --- sets particles boundary conditions
    warp.top.pboundxy = warp.periodic      #reflect#absorb
    warp.top.pbound0  = warp.absorb
    warp.top.pboundnz = warp.absorb

    if dim=="1d":
      warp.w3d.boundxy = warp.periodic
      warp.top.pboundxy = warp.periodic

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
    tmaxlab = (Lplasma_lab)/warp.clight
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

    #-------------------------------------------------------------------------------
    # set particles
    #-------------------------------------------------------------------------------
    weight     = 1.*dens0*warp.w3d.dx*warp.w3d.dy*warp.w3d.dz/(nppcellx*nppcelly*nppcellz) # weight of plasma macro-particles
    weightbeam = 0. # needs to be fixed

    # --- create plasma electron species
    elec = species.Species(type=species.Electron,weight=weight)
    # --- create plasma electron species
    prot = species.Species(type=species.Proton,weight=weight)
    # --- create e- beam species
    if l_beam:
      beam = species.Species(type=species.Electron,weight=weightbeam)
    # --- create plasma ion species
    # --- in this example, initial charge state is +5
    # --- charge states +6 and +7 are also considered
    ielec = species.Species(type=species.Electron,weight=weight)
    if l_ions:
      nions = 2
      ions=[]
      for i in range(nions):
        ions.append(species.Species(type=species.Krypton,charge_state=i+8,weight=weight,name='ions %g'%i))

    warp.top.pgroup.ndts=warp.nint(1./dtcoef)
    ielec.ndts=1

    #warp.top.pgroup.sm*=1000000

    #  ions.sm*=1.e20  # to give virtually infinite mass to ions
    warp.top.wpid = warp.nextpid() # creates data space for variable weights (needed for plasma ramps)
    warp.top.depos_order[...] = warp.top.depos_order[0,0] # sets deposition order of all species = those of species 0
    warp.top.efetch[...] = warp.top.efetch[0] # same for field gathering
    if dim in ["1d","2d"]:
      warp.top.depos_order[1,:]=1
    if dim=="1d":
      warp.top.depos_order[0,:]=1

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
    warp.top.fstype = -1 # sets field solver to None (desactivates electrostatic solver)
    warp.package('w3d')
    warp.generate()

    #-------------------------------------------------------------------------------
    # set a few shortcuts
    #-------------------------------------------------------------------------------
    pg = warp.top.pgroup

    #-------------------------------------------------------------------------------
    # sets tunnel ionizations
    #-------------------------------------------------------------------------------
    tunnel_ioniz = tunnel_ionization.TunnelIonization(stride=1)
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
    dx = warp.w3d.dx/nppcellx
    dy = warp.w3d.dy/nppcelly
    dz = warp.w3d.dz/nppcellz
    nx = nppcellx*warp.w3d.nx
    ny = nppcelly*warp.w3d.ny
    nz = nppcellz

    if dim in ["1d","2d"]:
      if dim=="1d":
        nx=1
        xp0,zp0 = warp.getmesh2d(0.,dx,0,
                            dz/2,dz,nz-1)
        yp0 = xp0*0.
      else:
        xp0,zp0 = warp.getmesh2d(warp.w3d.xmmin+dx/2,dx,nx-1,
                            -warp.w3d.dz+dz/2,dz,nz-1)
        yp0 = xp0*0.
    else:
      xp0,yp0,zp0 = warp.getmesh3d(warp.w3d.xmmin+dx/2,dx,nx-1,
                           warp.w3d.ymmin+dy/2,dy,ny-1,
                           -warp.w3d.dz+dz/2,dz,nz-1)

    zp0-=warp.minnd(zp0) # ensures that zp0 starts at 0

    # --- transform to 1D arrays
    xp0=xp0.flatten()
    yp0=yp0.flatten()
    zp0=zp0.flatten()

    # --- select particles within computational box of local processor
    ii=warp.compress((xp0>=warp.w3d.xmminlocal) & (xp0<warp.w3d.xmmaxlocal) & \
                (yp0>=warp.w3d.ymminlocal) & (yp0<warp.w3d.ymmaxlocal),numpy.arange(len(xp0)))
    xp0=warp.take(xp0,ii)
    yp0=warp.take(yp0,ii)
    zp0=warp.take(zp0,ii)

    # --- select particles within max_radius
    rp0=numpy.sqrt(xp0**2+yp0**2)
    ii=warp.compress(rp0<=max_radius,numpy.arange(len(rp0)))
    xp0=warp.take(xp0,ii)
    yp0=warp.take(yp0,ii)
    zp0=warp.take(zp0,ii)
    rp0=numpy.sqrt(xp0**2+yp0**2)

    # --- set the transverse profile
    wp0=plasma_trans_profile(rp0)

    if l_plasma:
        warp.installuserinjection(loadplasma)

    #-------------------------------------------------------------------------------
    # set laser pulse numpy.shape
    #-------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------
    # set laser amplitude by combining the pulse numpy.shape, laser profile, and laser phase
    #-------------------------------------------------------------------------------

    #------------------------------ add second laser--------------------------------#----------------set the time delay between two lasers--------------------------
    T_generate_laser=laser_risetime-laser_risetime1#2*laser_dura_fwhm  # time for generate laser pulses
    Tstart_laser1=T_generate_laser
    zstart_laser1=zstart_laser-21.25*lambda_laser
    vg=math.sqrt(1-kplab*kplab/k1lab/k1lab)*warp.clight

    #-------------------------------------------------------------------------------
    # initializes main field solver block
    #-------------------------------------------------------------------------------
    laser_func_dict={}
    laser_func_dict[0]=laser_func
    if not l_external_field:
      laser_func_dict[1]=laser_func1

    em = warp.EM3D(  laser_func=laser_func_dict,    # laser_func=laser_func,
                     laser_source_z=0.,
                     laser_source_v=-betafrm*warp.clight,
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

    if l_external_field:
        warp.installothereuser(add_external_laser)

    #-------------------------------------------------------------------------------
    # restarts from dump file
    #-------------------------------------------------------------------------------
    if l_restart:
        restore(dump_file)

    # --- load diagnostics
    #execfile('lpa_basic_diags.py')
    #execfile( diagnosticsScript )

    #-------------------------------------------------------------------------------
    # intializes e- beam
    #-------------------------------------------------------------------------------
    if l_beam:
      # --- add beam particles
      np_beam = 4000
      warp.top.vbeam = E_BEAM_BETA*warp.clight
      if me==0: # --- do only if processor 0
       if np_beam==1:
        # --- loads single test electron
        beam.addpart(x=array([0.]),
                     y=array([0.]),
                     z=array([0.]),
                     vx=array([0.]),
                     vy=array([0.]),
                     vz=array([E_BEAM_GAMMA*E_BEAM_BETA*warp.clight]),
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
           uz = ff.dp[::svstride]*GAMMAVZSIGMA+E_BEAM_GAMMA*E_BEAM_BETA*warp.clight
           gi = 1./math.sqrt(1.+(ux**2+uy**2+uz**2)/warp.clight**2)
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
                                  vzmean=E_BEAM_GAMMA*E_BEAM_BETA*warp.clight,
                                  lmomentum=True,
                                  lallindomain=True,
                                  zdist='regular')
      np_beam = beam.getn()
      # --- sets e- beam macro-particles weight
      if dim=="1d":
        beam.sw = (E_BEAM_DENSITY_PEAK*((E_BEAM_RADIUS*2)*warp.w3d.dx*warp.w3d.dy*E_BEAM_LENGTH*(2.*math.pi)))/np_beam
      if dim=="2d":
        beam.sw = (E_BEAM_DENSITY_PEAK*((E_BEAM_RADIUS*2)*warp.w3d.dy*E_BEAM_LENGTH*(2.*math.pi)))/np_beam
      if dim=="3d":
        beam.sw = (E_BEAM_DENSITY_PEAK*((E_BEAM_RADIUS*2)**2*E_BEAM_LENGTH*(2.*math.pi)*1.5))/np_beam
      # --- install e- beam diagnostic routine
      if sum(pg.nps)>0:
        # --- set beam position in lab frame
        pg.zp += zstart0-E_BEAM_DISTANCE_BEHIND_LASER-0.5*laser_total_length_lab
      # --- transform particle positions and velocity to boosted frame
      bf=Boosted_Frame(gammafrm,l_setselfb=0)
      zinit = getz().copy()
      bf.boost(beam,l_inject_plane=l_injectplane,lallindomain=1,l_rmzmean=0,zinject=-5.*warp.w3d.dz)
      particleboundaries3d(warp.top.pgroup,-1,False)

    #-------------------------------------------------------------------------------
    # register solver
    #-------------------------------------------------------------------------------
    print 'register solver'
    warp.registersolver(em)

    #-------------------------------------------------------------------------------
    # sets moving window velocity
    #-------------------------------------------------------------------------------
    if l_moving_window:
      #warp.top.vbeamfrm=BETAG_LINEAR*warp.clight
      warp.top.vbeamfrm=warp.clight
    #-------------------------------------------------------------------------------
    # set a few shortcuts
    #-------------------------------------------------------------------------------
    #el=elec
    #f = em.fields
    #bhist=numpy.ones([nions,warp.w3d.nz+1])
    #bhist[1:,:]=0.
    #def accuhist():
     # global bhist
     # for i in range(nions-2,-1,-1):
      #  Ex = em.gatherex()
      #  Ey = em.gatherey()
       # if me==0:
        #  E = math.sqrt(Ex*Ex+Ey*Ey)
         # E0 = numpy.where(E<=1.e-10*Eamp,1.e-10*Eamp,E)
         # dn = numpy.where(E==0.,0.,top.dt*tunnel_ioniz.GetADKrateSI(abs(E0),i+8,ions[i].type,top.dt))
         # bhist[i+1,:]+=dn*bhist[i]
         # bhist[i,:]-=dn*bhist[i]

    #installafterstep(accuhist)

    #window(1,hcp='lineout.cgm',dump=1,display='')
    #window(2,hcp='Ex.cgm',dump=1,display='')

    #if enableWarpLivePlots :
    #    installafterstep(liveplots)

    # make it visible to visit
    warp.listofallspecies = species.listofallspecies

    # add filtered species for ionized electrons
    # for now it should be at index 2...
    iElec = warp.listofallspecies[2]
    iElecCore = WarpVisItSpeciesFilters.CreateFilteredSpecies(iElec,'HaloFilter','halo')
    warp.listofallspecies.append(iElecCore)

    # print species
    i=0
    for s in warp.listofallspecies:
        print 'species %d %s %s'%(i,s.type.name,str(s.name))
        i += 1

    if (len(warp.listofallspecies) < 1):
        raise RuntimeError('no particle species!')

    return
