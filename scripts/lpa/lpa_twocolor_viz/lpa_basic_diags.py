from warp import *
import random
# --- define shortcut
ppg=ppgeneric

def labdata(z,t,ux,uy,uz,gi,uzfrm):
  if me==0 and l_verbose:print 'enter labdata'
  np = shape(z)[0]
  gammafrm = sqrt(1.+uzfrm**2/clight**2)
  zpr = gammafrm*z-uzfrm*t
  tpr = gammafrm*t-uzfrm*z/clight**2
  setu_in_uzboosted_frame3d(np,ux,uy,uz,gi,uzfrm,gammafrm)
  if me==0 and l_verbose:print 'exit labdata'
  return zpr,tpr,ux,uy,uz

def pzx(msize=1,color=red,titles=1,xscale=1.,yscale=1.,view=1):
  """
  plots ZX projection of e- beam
  """
  if not l_beam:return
  ppzx(msize=msize,color=color,titles=titles,xscale=xscale,yscale=yscale,view=view)
  try:
    ppzx(pgroup=bf.pgroup,color=blue,msize=msize,titles=titles,xscale=xscale,yscale=yscale,view=view)
  except:
    pass
  
def pzy(msize=1,color=red,titles=1,xscale=1.,yscale=1.,view=1):
  """
  plots ZY projection of e- beam
  """
  if not l_beam:return
  ppzy(msize=msize,color=color,titles=titles,xscale=xscale,yscale=yscale,view=view)
  try:
    ppzy(pgroup=bf.pgroup,color=blue,msize=msize,titles=titles,xscale=xscale,yscale=yscale,view=view)
  except:
    pass
  
def pxy(msize=1,color=red,titles=1,xscale=1.,yscale=1.,view=1):
  """
  plots XY projection of e- beam
  """
  if not l_beam:return
  ppxy(msize=msize,color=color,titles=titles,xscale=xscale,yscale=yscale,view=view)
  try:
    ppxy(pgroup=bf.pgroup,color=blue,msize=msize,titles=titles,xscale=xscale,yscale=yscale,view=view)
  except:
    pass
  
def pzxex(msize=1,view=1):
  em.pfex(l_transpose=dim<>'1d',direction=1,view=view)
  if not l_beam:return
  if dim=='1d':
    z=getz()
    if me==0:
      ppg(z*0,z,color=red,msize=msize,view=view)
    try:
      z=getz(pgroup=bf.pgroup)
      if me==0:
        ppg(z*0,z,color=blue,msize=msize,view=view)
    except:
      pass
  else:  
    pzx(msize=msize,titles=0,view=view)

def pzxey(msize=1,view=1):
  em.pfey(l_transpose=dim<>'1d',direction=1,view=view)
  if not l_beam:return
  if dim=='1d':
    z=getz()
    if me==0:
      ppg(z*0,z,color=red,msize=msize,view=view)
    try:
      z=getz(pgroup=bf.pgroup)
      if me==0:
        ppg(z*0,z,color=blue,msize=msize,view=view)
    except:
      pass
  else:  
    pzx(msize=msize,titles=0,view=view)

def pzxez(msize=1,view=1):
  em.pfez(l_transpose=dim<>'1d',direction=1,view=view)
  if not l_beam:return
  if dim=='1d':
    z=getz()
    if me==0:
      ppg(z*0,z,color=red,msize=msize,view=view)
    try:
      z=getz(pgroup=bf.pgroup)
      if me==0:
        ppg(z*0,z,color=blue,msize=msize,view=view)
    except:
      pass
  else:  
    pzx(msize=msize,titles=0,view=view)
  
def pzyex(msize=1,view=1):
  em.pfex(l_transpose=1,direction=0,view=view)
  pzy(msize=msize,titles=0,view=view)

def pzyey(msize=1,view=1):
  em.pfey(l_transpose=1,direction=0,view=view)
  pzy(msize=msize,titles=0,view=view)

def pzyez(msize=1,view=1):
  em.pfez(l_transpose=1,direction=0,view=view)
  pzy(msize=msize,titles=0,view=view)
  
def pxyex(msize=1,view=1):
  em.pfex(l_transpose=1,direction=2,view=view)
  pxy(msize=msize,titles=0,view=view)

def pxyey(msize=1,view=1):
  em.pfey(l_transpose=1,direction=2,view=view)
  pxy(msize=msize,titles=0,view=view)

def pxyez(msize=1,view=1):
  em.pfez(l_transpose=1,direction=2,view=view)
  pxy(msize=msize,titles=0,view=view)
  
zstart0lab=0.
dzstations=Lplasma_lab/nzstations
beamzstations = zstart0lab+arange(0.,Lplasma_lab,dzstations) # list of diag stations z-locations in lab frame

ekstations = zeros(shape(beamzstations),'d')
ppzstations = zeros(shape(beamzstations),'d')
xbarstations = zeros(shape(beamzstations),'d')
xpbarstations = zeros(shape(beamzstations),'d')
xsqstations = zeros(shape(beamzstations),'d')
xpsqstations = zeros(shape(beamzstations),'d')
xxpstations = zeros(shape(beamzstations),'d')
if dim == "3d":
  ybarstations = zeros(shape(beamzstations),'d')
  ypbarstations = zeros(shape(beamzstations),'d')
  ysqstations = zeros(shape(beamzstations),'d')
  ypsqstations = zeros(shape(beamzstations),'d')
  yypstations = zeros(shape(beamzstations),'d')
tbarstations = zeros(shape(beamzstations),'d')
tsqstations = zeros(shape(beamzstations),'d')
ekstationstime = zeros(shape(beamzstations),'d')
ekstationscnt = zeros(shape(beamzstations),'d')
ekstationscnt2 = zeros(shape(beamzstations),'d')
npz = 1001
pzbeamstations = (200.e6/dfact*arange(npz))/(npz-1)
pzstations = zeros((shape(beamzstations)[0],npz),'d')

top.zoldpid=nextpid()
def updatebeamstations():
  global timestart
#  if top.it%10<>0:return
  if me==0 and l_verbose:print 'enter updatebeamstations'
  # --- compute beta*gamma*c
  uzfrm=-betafrm*gammafrm*clight
  # --- get nb particles on each CPU
  np = getn(gather=0,bcast=0) 
  if np>0:
    # --- get z on each CPU
    z=getz(gather=0,bcast=0).copy()
    zold=getpid(id=top.zoldpid-1,gather=0,bcast=0)
    zoldlab = gammafrm*zold-uzfrm*(top.time-top.dt)
    # --- get z, time and velocities in lab frame
    zlab,tlab,uxlab,uylab,uzlab = labdata(z,
                                          top.time,
                                          getux(gather=0,bcast=0).copy(),
                                          getuy(gather=0,bcast=0).copy(),
                                          getuz(gather=0,bcast=0).copy(),
                                          getgaminv(gather=0,bcast=0).copy(),
                                          uzfrm=uzfrm)
    w = abs(zlab-zoldlab)/dzstations
    # --- get x,y on each CPU
    x = getx(gather=0,bcast=0).copy()
    y = gety(gather=0,bcast=0).copy()
    # --- compute gamma in lab frame
    myglab = sqrt(1.+(uxlab**2+uylab**2+uzlab**2)/clight**2)
    # --- compute kinetic energy in lab frame
    mykelab = beam.sm*(myglab-1.)*clight**2/echarge
    # --- defines cutoffs if particle selection is ON
    if l_pselect:
      # --- set threshold on radius
      XYcutoff = E_BEAM_RADIUS*5.
      # --- set threshold on longitudinal velocity
      UZcutoff = 0.95*E_BEAM_GAMMA*E_BEAM_BETA*clight
      # --- set threshold on energy
      KEcutoff = 1.e6 # eV
#      KEcutoff = None
#      XYcutoff = None
#      UZcutoff = None
    else:
      XYcutoff = None
      UZcutoff = None
      KEcutoff  = None
    if XYcutoff is not None:
      # --- select particle based on radius
      if dim=="3d":
        r2 = x*x+y*y
        XYcutoff2 = XYcutoff**2
        ii = compress(r2<XYcutoff2,arange(np))
      else:
        ii = compress(abs(x)<XYcutoff,arange(np))
      # --- get # of selected particles
      np = len(ii)
      # --- get weight, position, time, velocity and energy of selected particles
      w = take(w,ii)
      x = take(x,ii)
      y = take(y,ii)
      zlab = take(zlab,ii)
      tlab = take(tlab,ii)
      uxlab = take(uxlab,ii)
      uylab = take(uylab,ii)
      uzlab = take(uzlab,ii)
      mykelab = take(mykelab,ii)
    if UZcutoff is not None:
      # --- select particle based on longitudinal velocity
      ii = compress(uzlab>UZcutoff,arange(np))
      # --- get # of selected particles
      np = len(ii)
      # --- get weight, position, time, velocity and energy of selected particles
      w = take(w,ii)
      x = take(x,ii)
      y = take(y,ii)
      zlab = take(zlab,ii)
      tlab = take(tlab,ii)
      uxlab = take(uxlab,ii)
      uylab = take(uylab,ii)
      uzlab = take(uzlab,ii)
      mykelab = take(mykelab,ii)
    if KEcutoff is not None:
      # --- select particle based on energy
      ii = compress(mykelab>KEcutoff,arange(np))
      # --- get # of selected particles
      np = len(ii)
      # --- get weight, position, time, velocity and energy of selected particles
      w = take(w,ii)
      x = take(x,ii)
      y = take(y,ii)
      zlab = take(zlab,ii)
      tlab = take(tlab,ii)
      uxlab = take(uxlab,ii)
      uylab = take(uylab,ii)
      uzlab = take(uzlab,ii)
      mykelab = take(mykelab,ii)
    if np>0:
      xplab = uxlab/clight # normalized (gamma*beta*xp)
      yplab = uylab/clight # normalized (gamma*beta*yp) 
      nz = shape(ekstations)[0]
      deposgrid1dw(1,np,zlab,mykelab,w,nz-1,ekstations,ekstationscnt,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,beam.sm*uzlab*clight/echarge,w,nz-1,ppzstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,x,w,nz-1,xbarstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,x**2,w,nz-1,xsqstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,xplab,w,nz-1,xpbarstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,xplab**2,w,nz-1,xpsqstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,x*xplab,w,nz-1,xxpstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      if dim == "3d":
        deposgrid1dw(1,np,zlab,y,w,nz-1,ybarstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
        deposgrid1dw(1,np,zlab,y**2,w,nz-1,ysqstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
        deposgrid1dw(1,np,zlab,yplab,w,nz-1,ypbarstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
        deposgrid1dw(1,np,zlab,yplab**2,w,nz-1,ypsqstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
        deposgrid1dw(1,np,zlab,y*yplab,w,nz-1,yypstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,tlab,w,nz-1,tbarstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,tlab**2,w,nz-1,tsqstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      setgrid2dw(np,zlab,uzlab*beam.sm*clight/echarge,w,nz-1,npz-1,pzstations,
                  beamzstations[0],beamzstations[-1],pzbeamstations[0],pzbeamstations[-1])
  if top.it%hist_freq==0:
    savebeamstations()
  if me==0 and l_verbose:print 'exit updatebeamstations'

if l_beam:installafterstep(updatebeamstations)

def savebeamstations():
 if me==0 and l_verbose:print 'enter savebeamstations'
 pnums = parallelsum(ekstationscnt)
 pnum = where(pnums==0.,1.,pnums)
 ekst = parallelsum(ekstations)/pnum
 ppzst = parallelsum(ppzstations)/pnum
 xbar = parallelsum(xbarstations)/pnum
 xsq = parallelsum(xsqstations)/pnum
 xpbar = parallelsum(xpbarstations)/pnum
 xpsq = parallelsum(xpsqstations)/pnum
 xxp = parallelsum(xxpstations)/pnum
 if dim == "3d":
   ybar = parallelsum(ybarstations)/pnum
   ysq = parallelsum(ysqstations)/pnum
   ypbar = parallelsum(ypbarstations)/pnum
   ypsq = parallelsum(ypsqstations)/pnum
   yyp = parallelsum(yypstations)/pnum
 tbar = parallelsum(tbarstations)/pnum
 tsq = parallelsum(tsqstations)/pnum
 wti  = parallelsum(ekstationstime)/pnum
 pzst = parallelsum(pzstations)#*beam.sw
 if me==0:
  os.system('mv -f ebeamstations.pdb ebeamstationsold.pdb')
  f = PW.PW('ebeamstations.pdb')
  f.ekstations = ekst
  f.ppzstations = ppzst
  f.xbarstations = xbar
  f.xrmsstations = sqrt(xsq)
  f.xpbarstations = xpbar
  f.xprmsstations = sqrt(xpsq)
  f.xemitnstations = sqrt((xsq-xbar*xbar)*(xpsq-xpbar*xpbar)-(xxp-xbar*xpbar)**2)
  if dim == "3d":
    f.ybarstations = ybar
    f.yrmsstations = sqrt(ysq)
    f.ypbarstations = ypbar
    f.yprmsstations = sqrt(ypsq)
    f.yemitnstations = sqrt((ysq-ybar*ybar)*(ypsq-ypbar*ypbar)-(yyp-ybar*ypbar)**2)
  f.tbarstations = tbar
  f.trmsstations = sqrt(tsq-tbar*tbar)
  f.ekstationstime = wti
  f.pzstations = pzst
  f.beamzstations = beamzstations-zstart0lab
  f.pzbeamstations = pzbeamstations
  f.pnumstations = pnums
  f.nx = w3d.nx
  f.ny = w3d.ny
  f.nz = w3d.nz
  f.time = top.time
  f.dt = top.dt
  f.it = top.it
  f.stencil=stencil
  f.dim=dim
  f.close()
  os.system('rm -f ebeamstationsold.pdb')
 if me==0 and l_verbose:print 'exit savebeamstations'
  
def plke(view=1):
 global kelab,pxlab,pylab,pzlab,zhlab
 if me==0 and l_verbose:print 'enter plke'
 ekcnt = parallelsum(ekstationscnt)
 ekcnt = where(ekcnt==0.,1.,ekcnt)
 ekst = parallelsum(ekstations)/ekcnt
 if me==0:
    plsys(view)
    pla(ekst*1.e-6,beamzstations*1.e3,color=red)
    ptitles('Energy (MeV)','z (mm)','')
 if me==0 and l_verbose:print 'exit plke'


def maximum(value,default):
  try:
    return max(value)
  except ValueError:
    return default 

def minimum(value,default):
  try:
    return min(value)
  except ValueError:
    return default 

def pldenlineout_old(view=0,width=3.0):
  plsys(view)
  de = elec.get_density()/dens0
  dp = prot.get_density()/dens0
  di = ions[0].get_density()/dens0
  di1 = ions[1].get_density()/dens0
  wi_0 = ions[0].getweights()
  wi_1 = ions[1].getweights()
  if me==0:
   plsys(view)
   z = top.zgrid+w3d.zmmin+arange(w3d.nz+1)*w3d.dz
   if dim=="1d":
     dplineout=dp
     dilineout=di
     delineout=de
     ionized_number=sum(wi_1)*1.e-12
     total_number=sum(wi_0)*1.e-12
   else:
     dplineout=dp[w3d.nx/2,:]
     dilineout=di[w3d.nx/2,:]
     dilineout1=di1[w3d.nx/2,:]
     delineout=de[w3d.nx/2,:]
     ionized_number=sumnd(wi_1)*1.e-6
     ionized_number1=sum(dilineout1)*dens0*w3d.dx*w3d.dz*1.e-6
     total_number1=sum(dilineout)*dens0*w3d.dx*w3d.dz*1.e-6
   pla(dplineout,z*1.e3,color=blue,width=width)
   pla(dilineout*30,z*1.e3,color=red,width=width)
   pla(delineout,z*1.e3,color=black,width=width)
   ptitles('ionized # =%g; %g%%'%(ionized_number,100*ionized_number1/total_number1),'z (mm)','total # =%g'%total_number1)
   ylimits(0,3.0)
  
total_number1_max=0.
def pldenlineout(view=0,width=3.0):
  global total_number1_max
  plsys(view)
 # if me==0:print 'pldenlineout 1'
  if dim in ['1d','2d']:
    nxwide = 3
    xmin=-(nxwide+1)*w3d.dx/2
    xmax=(nxwide+1)*w3d.dx/2
#    nxwide=w3d.nx
#    xmin=w3d.xmmin
#    xmax=w3d.xmmax
    de = elec.get_density(nx=nxwide*nppcellx+1,xmin=xmin,xmax=xmax)/dens0
    de1= elec.get_density(nx=nxwide+1,xmin=xmin,xmax=xmax)/dens0
    die = ielec.get_density(nx=nxwide*nppcellx+1,xmin=xmin,xmax=xmax)/dens0
    die1= ielec.get_density(nx=nxwide+1,xmin=xmin,xmax=xmax)/dens0
    dp = prot.get_density(nx=nxwide*nppcellx+1,xmin=xmin,xmax=xmax)/dens0
    di = ions[0].get_density(nx=nxwide*nppcellx+1,xmin=xmin,xmax=xmax)/dens0
    di1 = ions[1].get_density(nx=nxwide*nppcellx+1,xmin=xmin,xmax=xmax)/dens0
  if dim=='3d':
    nxwide = 3
    xmin=-(nxwide+1)*w3d.dx/2
    xmax=(nxwide+1)*w3d.dx/2
    nywide = 3
    ymin=-(nywide+1)*w3d.dy/2
    ymax=(nywide+1)*w3d.dy/2
#    nxwide=w3d.nx
#    xmin=w3d.xmmin
#    xmax=w3d.xmmax
    de = elec.get_density(nx=nxwide*nppcellx+1,xmin=xmin,xmax=xmax,ny=nywide*nppcelly+1,ymin=ymin,ymax=ymax)/dens0
    de1= elec.get_density(nx=nxwide+1,xmin=xmin,xmax=xmax,ny=nywide+1,ymin=ymin,ymax=ymax)/dens0
    die = ielec.get_density(nx=nxwide*nppcellx+1,xmin=xmin,xmax=xmax,ny=nywide*nppcelly+1,ymin=ymin,ymax=ymax)/dens0
    die1= ielec.get_density(nx=nxwide+1,xmin=xmin,xmax=xmax,ny=nywide+1,ymin=ymin,ymax=ymax)/dens0
    dp = prot.get_density(nx=nxwide*nppcellx+1,xmin=xmin,xmax=xmax,ny=nywide*nppcelly+1,ymin=ymin,ymax=ymax)/dens0
    di = ions[0].get_density(nx=nxwide*nppcellx+1,xmin=xmin,xmax=xmax,ny=nywide*nppcelly+1,ymin=ymin,ymax=ymax)/dens0
    di1 = ions[1].get_density(nx=nxwide*nppcellx+1,xmin=xmin,xmax=xmax,ny=nywide*nppcelly+1,ymin=ymin,ymax=ymax)/dens0
  wi_0 = ions[0].getweights(gather=0)
  wi_1 = ions[1].getweights(gather=0)
  if dim=="1d":
     ionized_number=ionized_number1=globalsum(wi_1)*1.e-12
     total_number=total_number1=globalsum(wi_0)*1.e-12
  if dim=="2d":
    xi_0 = ions[0].getx(gather=0)
    xi_1 = ions[1].getx(gather=0)
    ii = compress( (xi_0<0.5*w3d.dx/nppcellx) & (xi_0>-0.5*w3d.dx/nppcellx),arange(shape(xi_0)[0]))
    wi_0_axis = take(wi_0,ii)
    ii = compress( (xi_1<0.5*w3d.dx/nppcellx) & (xi_1>-0.5*w3d.dx/nppcellx),arange(shape(xi_1)[0]))
    wi_1_axis = take(wi_1,ii)
    ionized_number1 = globalsum(wi_1_axis)
    total_number1 = globalsum(wi_0_axis)
    ionized_number=globalsum(wi_1)*1.e-6
  if dim=="3d":
    xi_0 = ions[0].getx(gather=0)
    xi_1 = ions[1].getx(gather=0)
    yi_0 = ions[0].gety(gather=0)
    yi_1 = ions[1].gety(gather=0)
    ii = compress( (xi_0<0.5*w3d.dx/nppcellx) \
                 & (xi_0>-0.5*w3d.dx/nppcellx) \
                 & (yi_0<0.5*w3d.dy/nppcelly) \
                 & (yi_0>-0.5*w3d.dy/nppcelly) \
                 ,arange(shape(xi_0)[0]))
    wi_0_axis = take(wi_0,ii)
    ii = compress( (xi_1<0.5*w3d.dx/nppcellx) \
                 & (xi_1>-0.5*w3d.dx/nppcellx) \
                 & (yi_1<0.5*w3d.dy/nppcelly) \
                 & (yi_1>-0.5*w3d.dy/nppcelly) \
                 ,arange(shape(xi_1)[0]))
    wi_1_axis = take(wi_1,ii)
    ionized_number1 = globalsum(wi_1_axis)
    total_number1 = globalsum(wi_0_axis)
    ionized_number=globalsum(wi_1)*1.e-6
  if me==0:
   total_number1_max = max(total_number1,total_number1_max)
   plsys(view)
   z = top.zgrid+w3d.zmmin+arange(w3d.nz+1)*w3d.dz
   if dim=="1d":
     dplineout=dp
     dilineout=di
     delineout1=de1
     dielineout1=die1
   if dim=="2d":
     dplineout=dp[(nxwide*nppcellx+1)/2,:]
     dilineout=di[(nxwide*nppcellx+1)/2,:]
     dilineout1=di1[(nxwide*nppcellx+1)/2,:]
     delineout=de[(nxwide*nppcellx+1)/2,:]
     delineout1=de1[nxwide/2,:]
     dielineout=die[(nxwide*nppcellx+1)/2,:]
     dielineout1=die1[nxwide/2,:]
   if dim=="3d":
     dplineout=dp[(nxwide*nppcellx+1)/2,(nywide*nppcelly+1)/2,:]
     dilineout=di[(nxwide*nppcellx+1)/2,(nywide*nppcelly+1)/2,:]
     dilineout1=di1[(nxwide*nppcellx+1)/2,(nywide*nppcelly+1)/2,:]
     delineout=de[(nxwide*nppcellx+1)/2,(nywide*nppcelly+1)/2,:]
     delineout1=de1[nxwide/2,nywide/2,:]
     dielineout=die[(nxwide*nppcellx+1)/2,(nywide*nppcelly+1)/2,:]
     dielineout1=die1[nxwide/2,nywide/2,:]
#     ionized_number1=sum(dilineout1)*dens0*w3d.dx*w3d.dz*1.e-6
#     total_number1=sum(dilineout)*dens0*w3d.dx*w3d.dz*1.e-6
   pla(dplineout,z*1.e3,color=blue,width=width)
   pla(dilineout*2,z*1.e3,color=red,width=width)
   pla(delineout1,z*1.e3,color=black,width=width)
   pla(dielineout1,z*1.e3,color=green,width=width)
   if total_number1_max>0.:
     fraction = ionized_number1/total_number1_max
   else:
     fraction = 0.
   ptitles('ionized # =%g; %g; %g%%'%(ionized_number,ionized_number1*1.e-6,100*fraction),'z (mm)','')
  #ylimits(0,3.0)
  

#def pledensall(view=0,width=5):
 # plsys(view)
  #de = elec.get_density()*elec.charge/echarge/dens0
  #dp = prot.get_density()*prot.charge/echarge/dens0
  #di = ions[0].get_density()*ions[0].charge/echarge/dens0
  #z = top.zgrid+w3d.zmmin+arange(w3d.nz+1)*w3d.dz
  #pla(dp+di+de,z,color=black,width=width)
  #pla(dp,z,color=blue,width=width)
  #pla(di,z,color=green,width=width)
  #pla(de,z,color=red,width=width)
  #ptitles('charge density/dens0','z [m]','')

  
#t_his=AppendableArray()
#ex_his=AppendableArray(unitshape=shape(em.gatherex()))
#ey_his=AppendableArray(unitshape=shape(em.gatherey()))
#ez_his=AppendableArray(unitshape=shape(em.gatherez()))
#ey_his_freq=600
#def save_field_data():
  #ex=em.gatherex()
  #ey=em.gatherey()
  #ez=em.gatherez()
  #if me==0:
    #if top.it%ey_his_freq==0:
    # ex_his.append(ex)   
    # ey_his.append(ey)   
    # ez_his.append(ez)
    # t_his.append(top.time)   
    # z1=w3d.zmmin+arange(w3d.nz+1)*w3d.dz
    # f=PW.PW('field_data')
    # f.z=z1
    # f.t=t_his[...]
    # f.ex=ex_his[...]
    # f.ey=ey_his[...]
    # f.ez=ez_his[...]
    # f.close()


def plexezlineout(view=0,width=2.8):
  global em
  #ex=em.gatherex()
  #ez=em.gatherez()
  
  if me==0:
   plsys(view)
   z=top.zgrid+w3d.zmmin+arange(w3d.nz+1)*w3d.dz
   exlineout=ex[w3d.nx/2,:]
   ezlineout=ez[w3d.nx/2,:]
   pla(exlineout*1.e-9,z*1.e3,color=black,width=width)
   pla(ezlineout*1.e-9*20,z*1.e3,color=red,width=3.0)
   ptitles('Ex, Ez*20 (GV/m)','z (mm)','')
   #ylimits(-1300,1300)

#def plex(view=0):
 # em.pfex(l_transpose=1,xscale=1.e3,yscale=1.e3,cmin=-1.3e12,cmax=1.3e12,titles=0,view=view)
 # palette('bluewhitered.gp')
 # ptitles('Ex (V/m)','z (mm)','x (mm)')  

nzplots= min(w3d.nz,1000)
nxplots= min(w3d.nx,100)
nyplots= min(w3d.ny,100)
dxplots = (w3d.xmmax-w3d.xmmin)/nxplots
dyplots = (w3d.ymmax-w3d.ymmin)/nyplots
dzplots = (w3d.zmmax-w3d.zmmin)/nzplots
if dim=="1d":
    nxplots=1
    xplots,zplots = getmesh2d(0.,dxplots,0,
                              dzplots/2,dzplots,nzplots-1)
    yplots = xplots*0.
if dim in ["2d","3d"]:
    xplots,zplots = getmesh2d(w3d.xmmin+dxplots/2,dxplots,nxplots-1,
                              w3d.zmmin+dzplots/2,dzplots,nzplots-1)
    yplots = xplots*0.
#else:
#  xplots,yplots,zplots = getmesh3d(w3d.xmmin+dxplots/2,dxplots,nxplots-1,
#                          w3d.ymmin+dyplots/2,dyplots,nyplots-1,
#                          w3d.zmmin+dzplots/2,dzplots,nzplots-1)

# --- transform to 1D arrays
xplotsfl=xplots.flatten()
yplotsfl=yplots.flatten()
zplotsfl=zplots.flatten()

def plex2d(view=0):
  global em,xplots,yplots,zplots,xplotsfl,yplotsfl,zplotsfl
  plsys(view)
  ex,ey,ez,bx,by,bz=em.getfieldsfrompositions(xplotsfl,yplotsfl,zplotsfl+top.zgrid)
  if l_external_field:
    t = top.time-(zplotsfl+top.zgrid-zstart_laser1)/clight
    ex1,ey1=laser_func1(xplotsfl,yplotsfl,t)
    ex+=ex1
    ey+=ey1 
  ex=ex.reshape(shape(xplots))
  palette('field.gp')
  ppg(transpose(ex*1.e-12), \
      xmin=1.e3*(zplots[0,0]+top.zgrid-0.5*dzplots),\
      xmax=1.e3*(zplots[0,-1]+top.zgrid-0.5*dzplots),\
      ymin=1.e3*(xplots[0,0]-0.5*dxplots),\
      ymax=1.e3*(xplots[-1,0]-0.5*dxplots),\
      titles=0,view=view)
  ptitles('Ex (TV/m)','z (mm)','x (mm)')


def plex2dwake(view=0):
  global em,xplots,yplots,zplots,xplotsfl,yplotsfl,zplotsfl
  plsys(view)
  ex,ey,ez,bx,by,bz=em.getfieldsfrompositions(xplotsfl,yplotsfl,zplotsfl+top.zgrid)
  if l_external_field:
    t = top.time-(zplotsfl+top.zgrid-zstart_laser1)/clight
    ex1,ey1=laser_func1(xplotsfl,yplotsfl,t)
    ex+=ex1
    ey+=ey1
    by+=ex1/clight
  ex=ex.reshape(shape(xplots))
  by=by.reshape(shape(xplots))
  ez=ez.reshape(shape(xplots))
  force=ex-by*clight
  #palette('bluewhitered.gp')
  force_max=maximum(force,0)*1.e-9
  ppg(transpose(force*1.e-9), \
      xmin=1.e3*(zplots[0,0]+top.zgrid-0.5*dzplots),\
      xmax=1.e3*(zplots[0,-1]+top.zgrid-0.5*dzplots),\
      ymin=1.e3*(xplots[0,0]-0.5*dxplots),\
      ymax=1.e3*(xplots[-1,0]-0.5*dxplots),\
      titles=0,view=view)
  ptitles('Ex-c*By(GV/m)','z (mm)','x (mm)')
  field_his_freq=30800
  z1=top.zgrid+w3d.zmmin+arange(w3d.nz+1)*w3d.dz
  zmax=max(z1)*1.e3
  if me==0:
    if top.it%field_his_freq==0:
     f=PW.PW('field_data_z=%.2f'%zmax)
     f.ex=ex
     f.by=by
     f.ez=ez
     f.close() 
  



def plex(view=0,width=2.8):
  global em
#  ex=em.gatherex()
#  ey=em.gatherey()
#  ez=em.gatherez()
  nz = 10000
  zmax = w3d.zmmax+top.zgrid
  zmin = w3d.zmmin+top.zgrid
  dz = (zmax-zmin)/nz
  z = zmin+arange(nz)*dz
  x=y=z*0.
 # if me==0:print 'plex 1'
  ex,ey,ez,bx,by,bz=em.getfieldsfrompositions(x,y,z)
  if l_external_field:
    t = top.time-(z-zstart_laser1)/clight
    ex1,ey1=laser_func1(x,y,t)
    ex+=ex1
    ey+=ey1
#  if me==0:print 'plex 2'
  if me==0:
   plsys(view)
   ezmax=maximum(ez,0)*1.e-9
#   z=top.zgrid+w3d.zmmin+arange(w3d.nz+1)*w3d.dz
   pla(ex*1.e-9,z*1.e3,color=black,width=width)
   pla(ez*1.e-9*20,z*1.e3,color=red,width=3.0)
   ptitles('Maximum Ez(GV/m)=%g'%ezmax,'z (mm)','Ez*20')
   #ylimits(-1300,1300)
 # if me==0:print 'plex 3'




def plez(view=0):
  global em,xplots,yplots,zplots,xplotsfl,yplotsfl,zplotsfl
  plsys(view)
  ex,ey,ez,bx,by,bz=em.getfieldsfrompositions(xplotsfl,yplotsfl,zplotsfl+top.zgrid)
  #if l_external_field:
  #  t = top.time-(zplotsfl+top.zgrid-zstart_laser1)/clight
  #  ex1,ey1=laser_func1(xplotsfl,yplotsfl,t)
  #  ex+=ex1
  #  ey+=ey1 
  ez=ez.reshape(shape(xplots))
  ppg(transpose(ez*1.e-9), \
      xmin=1.e3*(zplots[0,0]+top.zgrid-0.5*dzplots),\
      xmax=1.e3*(zplots[0,-1]+top.zgrid-0.5*dzplots),\
      ymin=1.e3*(xplots[0,0]-0.5*dxplots),\
      ymax=1.e3*(xplots[-1,0]-0.5*dxplots),\
      cmin=-35,cmax=35,titles=0,view=view)
  ezmax = globalmax(em.getez())*1.e-9
  ptitles('Maximum Ez(GV/m)=%g'%ezmax,'z (mm)','x (mm)')

def pldens(view=0):
  plsys(view)
  de = elec.get_density(nx=min([w3d.nx,100]),nz=min([w3d.nz,200]))
  if me==0:
    de1=de.transpose()
  if ielec.getn()>0:
    die = ielec.get_density(nx=min([w3d.nx,100]),nz=min([w3d.nz,200]))
    if me==0:
      de1+=die.transpose()
  if me==0:
    ppg(de1, \
        xmin=w3d.zmmin+top.zgrid-0.5*w3d.dz, \
        xmax=w3d.zmmax+top.zgrid-0.5*w3d.dz, \
        ymin=w3d.xmmin-0.5*w3d.dx, \
        ymax=w3d.xmmax-0.5*w3d.dx, \
        xscale=1.e3,yscale=1.e3)
    palette('density.gp')
    ptitles('Electron density','z (mm)','x (mm)')

def pldens1d(view=0,width=3.0):
  plsys(view)
  de = elec.get_density()/dens0
  dp = prot.get_density()/dens0
  di = ions[0].get_density()/dens0
  wi_0 = ions[0].getweights()
  wi_1 = ions[1].getweights()
  ionized_number=sum(wi_1)*1.e-12
  total_number=sum(wi_0)*1.e-12
  if me==0:
   plsys(view)
   z = top.zgrid+w3d.zmmin+arange(w3d.nz+1)*w3d.dz
   pla(dp,z*1.e3,color=blue,width=width)
   pla(di*30,z*1.e3,color=red,width=width)
   pla(de,z*1.e3,color=black,width=width)
   ptitles('ionized # =%g'%ionized_number,'z (mm)','total # =%g'%total_number)
   #ylimits(0,3.0)

energy_trap=None  
energy_trapy=None  
n_trap=0
def plzuz(view=0,msize=1.5,nplot=100000):
  global energy_trap,energy,ux_trap,uy_trap,uz_trap,z_trap,emitt_trap,emitt_trapy,we_trap,x_trap,Ex_trap,rms_x2,rms_ux2,rms_xux2,By_trap,vz_trap,Ez_trap,z_plasma,x_plasma,Ex_plasma,By_plasma,vz_plasma,n_trap

  z1=top.zgrid+w3d.zmmin+arange(w3d.nz+1)*w3d.dz
  zmax=max(z1)
  zmin=min(z1)
  z2=zmax-(zmax-zmin)*9.0/14.0

  cond = dim <>"3d"
  if cond:
    Ex=elec.getex(gather=0,bcast=0)
    Ey=elec.getey(gather=0,bcast=0)
    Ez=elec.getez(gather=0,bcast=0)
    Bx=elec.getbx(gather=0,bcast=0)
    By=elec.getby(gather=0,bcast=0)
    Bz=elec.getbz(gather=0,bcast=0)
    vz=elec.getvz(gather=0,bcast=0)
    z=elec.getz(gather=0,bcast=0)
    x=elec.getx(gather=0,bcast=0)
    y=elec.gety(gather=0,bcast=0)
    ux=elec.getux(gather=0,bcast=0)/clight
    uy=elec.getuy(gather=0,bcast=0)/clight
    uz=elec.getuz(gather=0,bcast=0)/clight
    we=elec.getweights(gather=0,bcast=0)
    gamma=sqrt(1+ux*ux+uy*uy+uz*uz)
    betaz=uz/gamma
    betag=sqrt(1-wplab*wplab/w0lab/w0lab)  ####lab frame laser group velocity
    energy=(gamma-1)*0.511
    nlocal=elec.getn(gather=0,bcast=0) 
    n=globalsum(nlocal)
    if n==0:return
    # --- get randomly a max of around nplot particles to plot
    ii = arange(nlocal)
    frac = min(1.,float(nplot)/n)
    random.shuffle(ii)
    ntoplot = int(frac*nlocal)
    ii = ii[:ntoplot]  
    zplot = gatherarray(take(z,ii))
    uzplot = gatherarray(take(uz,ii))
    if me==0:
      plsys(view)
    #ppgeneric(uzplot,z*1.e3,zzplot=we_trap,color='density',nx=20,ny=20,contours=100,ncolor=100,filled=1,lframe=1,pplimits=(zmin*1.e3,zmax*1.e3,'e','e'))
      ppgeneric(uzplot,zplot*1.e3,color=black,msize=msize,lframe=1,pplimits=(zmin*1.e3,zmax*1.e3,'e','e'))
#ptitles('Charge [pC/(um*um)]=%g'%charge_trap,'z (mm)','uz')
  # --- trap electrons
  #ii=compress((betaz>0) & (betaz>=betag) & (z<z2),arange(n))
  Ex=ielec.getex(gather=0,bcast=0)
  Ey=ielec.getey(gather=0,bcast=0)
  Ez=ielec.getez(gather=0,bcast=0)
  Bx=ielec.getbx(gather=0,bcast=0)
  By=ielec.getby(gather=0,bcast=0)
  Bz=ielec.getbz(gather=0,bcast=0)
  vz=ielec.getvz(gather=0,bcast=0)
  z=ielec.getz(gather=0,bcast=0)
  x=ielec.getx(gather=0,bcast=0)
  y=ielec.gety(gather=0,bcast=0)
  ux=ielec.getux(gather=0,bcast=0)/clight
  uy=ielec.getuy(gather=0,bcast=0)/clight
  uz=ielec.getuz(gather=0,bcast=0)/clight
  we=ielec.getweights(gather=0,bcast=0)
  nlocal=ielec.getn(gather=0,bcast=0) 
  if nplot==0:return
  zplot = gatherarray(z)
  uzplot = gatherarray(uz)
  if cond:
    if me==0:
      ppgeneric(uzplot,zplot*1.e3,color=red,msize=msize,lframe=1,pplimits=(zmin*1.e3,zmax*1.e3,'e','e'))
  gamma=sqrt(1+ux*ux+uy*uy+uz*uz)
  betaz=uz/gamma
  betag=sqrt(1-wplab*wplab/w0lab/w0lab)  ####lab frame laser group velocity
  energy=(gamma-1)*0.511

  ii=compress(uz>=0.8,arange(nlocal))
  n_trap_local=len(ii)
  n_trap = globalsum(n_trap_local)
  if n_trap==0:return
  ux_trap=gatherarray(take(ux,ii))
  uy_trap=gatherarray(take(uy,ii))
  uz_trap=gatherarray(take(uz,ii))
  z_trap=gatherarray(take(z,ii))
  x_trap=gatherarray(take(x,ii))
  y_trap=gatherarray(take(y,ii))
  Ex_trap=gatherarray(take(Ex,ii))
  Ey_trap=gatherarray(take(Ey,ii))
  Ez_trap=gatherarray(take(Ez,ii))
  Bx_trap=gatherarray(take(Bx,ii))
  By_trap=gatherarray(take(By,ii))
  Bz_trap=gatherarray(take(Bz,ii))
  vz_trap=gatherarray(take(vz,ii)) 
  we_trap=gatherarray(take(we,ii))
  if me>0:return
  xvar = sqrt(var(x_trap))
  zvar = sqrt(var(z_trap))
  uxvar = sqrt(var(ux_trap))
  uzvar = sqrt(var(uz_trap))
  xave=ave(x_trap)
  zave=ave(z_trap)
  uxave=ave(ux_trap)
  uzave=ave(uz_trap)
  gamma_trap=sqrt(1+ux_trap*ux_trap+uy_trap*uy_trap+uz_trap*uz_trap)
  energy_trap=(gamma_trap-1)*0.511
  #charge_trap=we_trap*echarge
  #charge_trap=sum(charge_trap)
  number_trap=sum(we_trap)*1.e-6
#  ztrapmax=max(z_trap)
#  ztrapmin=min(z_trap)
#  xtrapmax=max(x_trap)
#  xtrapmin=min(x_trap)
#  ij=compress((betaz<betag) & (z>=ztrapmin) & (z<=ztrapmax) & (x>=xtrapmin) & (x<=xtrapmax),arange(n))
 # ij=compress(z<z2,arange(n))
 # z_plasma=gatherarray(take(z,ij))
 # x_plasma=gatherarray(take(x,ij))
 # Ex_plasma=gatherarray(take(Ex,ij))
 # By_plasma=gatherarray(take(By,ij))
 # Ez_plasma=gatherarray(take(Ez,ij))
 # ux_plasma=gatherarray(take(ux,ij))
 # uy_plasma=gatherarray(take(uy,ij))
 # uz_plasma=gatherarray(take(uz,ij))
 # we_plasma=gatherarray(take(we,ij))
  rms_x2=sum(x_trap*x_trap*we_trap)/sum(we_trap)
  rms_ux2=sum(ux_trap*ux_trap*we_trap)/sum(we_trap)
  rms_xux=sum(x_trap*ux_trap*we_trap)/sum(we_trap)
  rms_xux2=rms_xux*rms_xux
  emitt_trap=sqrt(rms_x2*rms_ux2-rms_xux2)
  rms_y2=sum(y_trap*y_trap*we_trap)/sum(we_trap)
  rms_uy2=sum(uy_trap*uy_trap*we_trap)/sum(we_trap)
  rms_yuy=sum(y_trap*uy_trap*we_trap)/sum(we_trap)
  rms_yuy2=rms_yuy*rms_yuy
  emitt_trapy=sqrt(rms_y2*rms_uy2-rms_yuy2)
  ele_his_freq=1#30800
  zmax1=zmax*1.e3
  if me==0:
    if top.it%ele_his_freq==0:
     f=PW.PW('electron_data_z=%.2f'%zmax1)
     f.xtrap=x_trap
     f.ytrap=y_trap
     f.ztrap=z_trap
     f.uxtrap=ux_trap
     f.uytrap=uy_trap
     f.uztrap=uz_trap
     f.wetrap=we_trap
     f.Extrap=Ex_trap
     f.Eytrap=Ey_trap
     f.Eztrap=Ez_trap
     f.Bxtrap=Bx_trap
     f.Bytrap=By_trap
     f.Bztrap=Bz_trap
    f.close() 
  if me==0:
    ptitles('Trapped # [/um]=%g'%number_trap,'z (mm)','uz')

def plspectrum(view=0):
   global energy_trap, we_trap, n_trap
   if n_trap==0:return
   if me==0:
    we_trap1=we_trap*1.e-6
    he,be=histogram(energy_trap,20,weights=we_trap1)
    plsys(view)
    energymax=maximum(energy_trap,0)
    rms_E2=sum(energy_trap*energy_trap*we_trap)/sum(we_trap)
    rms_E=sum(energy_trap*we_trap)/sum(we_trap)
    rms_E22=rms_E*rms_E
    rms_deltaE=sqrt(rms_E2-rms_E22)
    Energy_spread_rms=rms_deltaE/rms_E
    plh(he,be)
    ptitles('%g; %g; %g%%'%(energymax,rms_deltaE,Energy_spread_rms*100),'Energy (MeV)','')


def plxux(view=0,msize=4.0):
  global energy_trap,ux_trap,x_trap,z_trap,we_trap,emitt_trap,rms_x2,rms_ux2,rms_xux2
  if energy_trap is None:return
  if me==0:
   plsys(view)
 #  z1=top.zgrid+w3d.zmmin+arange(w3d.nz+1)*w3d.dz
 #  zmax=max(z1)
 #  zmin=min(z1)+4.e-5
 #  zcenter=(zmax+zmin)/2.0-9.e-5
   ux_max=maximum(ux_trap,0)
   #ppgeneric(ux_trap,x_trap*1.e6,color=black,msize=msize)
   palette('field.gp')
   ppgeneric(ux_trap,x_trap*1.e6,weights=we_trap*1.e-6,msize=4.0,color='density',nx=50,ny=50,ncolor=200)
   ptitles('%.4f; %g; %g; %g'%(emitt_trap*1.e6,sqrt(rms_x2)*1.e6,sqrt(rms_ux2),sqrt(rms_xux2)*1.e6),'x (um)','')  

def plex_x_z(view=0,msize=4):
  global energy_trap,Ex_trap,x_trap,z_trap,we_trap,By_trap,vz_trap
  if energy_trap is None:return
  if me==0:
   plsys(view)
   Ex_max=maximum(Ex_trap,0)*1.e-9
   force_trap=vz_trap*By_trap-Ex_trap
   ppg(force_trap*1.e-9,x_trap*1.e6,weights=z_trap*1.e6,msize=4,ncolor=200,color='density')
   ptitles('Max Ex(GV/m)=%g'%Ex_max,'x (um)','-(Ex-Vz*By) (GV/m)')


def plx_z_ex(view=0,msize=4):
  global energy_trap,Ex_trap,x_trap,z_trap,we_trap,By_trap,vz_trap
  if energy_trap is None:return
  if me==0:
   plsys(view)
   rms_z2=sum(z_trap*z_trap*we_trap)/sum(we_trap)
   rms_z=sqrt(rms_z2)*1.e6
   force_trap=vz_trap*By_trap-Ex_trap
   ppg(x_trap*1.e6,z_trap*1.e6,weights=force_trap*1.e-9,msize=4,ncolor=200,color='density')
   ptitles('%g'%rms_z,'z (um)','x (um)')



#def plx_z_ex_p(view=0,msize=4):
#  global energy_trap,Ex_plasma,x_plasma,z_plasma,By_plasma,vz_plasma
#  if energy_trap is None:return
#  if me==0:
#   plsys(view)
#   force_plasma=vz_plasma*By_plasma-Ex_plasma
#   ppg(x_plasma*1.e6,z_plasma*1.e6,weights=force_plasma*1.e-9,msize=4,ncolor=200,color='density')
#   ptitles('','z (um)','x (um)')


z_propag=[]
emitx_hist=[]
emity_hist=[]
def save_emittance():
   global emitt_trap,emitt_trapy,energy_trap, z_trap, z_propag, emitx_hist, emity_hist
   if energy_trap is None:return
   if me==0:
#      z_propag.append(top.time*clight*1.e3)
      z_propag.append(ave(z_trap)*1.e3)
      emitx_hist.append(emitt_trap)
      emity_hist.append(emitt_trapy)
      #pla(emitt_trap,z_propag)
      #ppg(emitt_trap,z_propag,msize=5.0)
      f=PW.PW('emittance_data')
      f.z=z_propag
      f.emitx=emitx_hist
      f.emity=emity_hist
      f.close()


#def plzuy(view=0,msize=1.5):
 # if me==0:
 #  plsys(view)
  # z1=top.zgrid+w3d.zmmin+arange(w3d.nz+1)*w3d.dz
  # zmax=max(z1)
  # zmin=min(z1)+4.e-5
  # zcenter=(zmax+zmin)/2.0-9.e-5
  # uy_max=maximum(uy_trap,0)
  # ppgeneric(uy_trap,z_trap*1.e3,color=black,msize=msize,lframe=1,pplimits=(zmin*1.e3,zcenter*1.e3,'e','e'))
  # ptitles('Maximum uy=%g'%uy_max,'z (mm)','uy')  

if me==0:
  tottime = AppendableArray()
  def accuttime():
    global tottime
    tottime.append(time.clock())
    if top.it%200==0:
      f=PW.PW('tottime.pdb')
      f.time=tottime[:]
      f.close()
  installafterstep(accuttime)

  
