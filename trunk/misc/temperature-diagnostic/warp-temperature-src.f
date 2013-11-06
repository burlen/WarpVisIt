c=============================================================================
      subroutine reset_temperature(is)
      use Temperatures
      use Timers
      integer(ISZ):: is
      real(kind=8):: timetemp,wtime
      timetemp = wtime()
        pnumt = 0.
        pnumtw = 0.
        vxbart = 0.
        vybart = 0.
        vzbart = 0.
        vxsqbart = 0.
        vysqbart = 0.
        vzsqbart = 0.
        kebart = 0.
        kesqbart = 0.
        xkebart = 0.
        ykebart = 0.
        zkebart = 0.
        tempxz(:,is) = 0.
        tempyz(:,is) = 0.
        tempzz(:,is) = 0.
        tempx(:,:,:,is) = 0.
        tempy(:,:,:,is) = 0.
        tempz(:,:,:,is) = 0.
        dke(:,:,:,is) = 0.
        if(l_temp_rmcorrelations .or. l_temp_rmcrosscorrelations) then
          xbart = 0.
          ybart = 0.
          zbart = 0.
          xsqbart = 0.
          ysqbart = 0.
          zsqbart = 0.
        endif
        if(l_temp_rmcorrelations) then
          xvxbart = 0.
          yvybart = 0.
          zvzbart = 0.
        endif
        if(l_temp_rmcrosscorrelations) then
          xybart = 0.
          xzbart = 0.
          yzbart = 0.
          xvybart = 0.
          xvzbart = 0.
          yvxbart = 0.
          yvzbart = 0.
          zvxbart = 0.
          zvybart = 0.
        endif
      temperaturestime = temperaturestime + (wtime() - timetemp)
      return
      end
c=============================================================================
      subroutine accumulate_temperature(np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,dt,
     &                          uxpo,uypo,uzpo,is,wp,lw,lrtheta,l2symtry,l4symtry)
c     Compute temperature in Z-slices for species 'is' on a 3-D grid:
c       - the slices can have any position and thickness but cannot overlap,
c       - the min and max of each slice in x, y and z are given respectively in the arrays
c         tslicexmin, tslicexmax, tsliceymin, tsliceymax, tslicezmin and tslicezmax
c         (the reason for having slices with different thickness and dimensions is
c          to allow the temperature measurement to the shape of a distribution, like
c          for example a beam extending over several quadrupoles and accelerating gaps),
c       - the x, y and z temperatures are given in the arrays tempx, tempy and tempz,
c         while averages in each slice are given in the arrays tempxz, tempyz, tempzz.
c       - the calculation is done in three parts:
c         o reset_temperature: zero out all moments
c         o accumulate_temperature: accumulate moments from particles
c         o finalize_temperature: compute final quantities
c       - l_temp_collapseinz=.true.: collapse slices in z, i.e. align Z-locations
c                                    of particles using current velocity (uxp,uyp,yzp) and
c                                    velocity from previous time step (uxpo,uypo,uzpo)
c       - lrtheta=.true.: radial and azimuthal are computed in place of X and Y,
c       - the default temperature unit is in electron-volt. To select the units,
c         set the variable t_unit to the default integers evolt, joule or kelvin.

      use Beam_acc
      use InDiag
      use Picglb
      use ExtPart
      use Temperatures
      use Timers

      integer(ISZ):: np,is
      real(kind=8):: w,dt,ke
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uxp(np), uyp(np), uzp(np), gaminv(np), wp(np)
      real(kind=8):: uxpo(np), uypo(np), uzpo(np) ! uxp, uyp anx uzp of previous time step
      logical(ISZ):: lw,lrtheta,l2symtry,l4symtry

      integer(ISZ):: ip,its,izl,ixt,iyt,izt
      real(kind=8):: dti,wt,ddx,ddy,oddx,oddy,xt,yt,wpp,z_local
      real(kind=8):: oneondt,clighti,vzi,zc
      real(kind=8):: xpt,ypt,zpt,vxpt,vypt,vzpt,cost,sint,rpt
      real(kind=8):: timetemp,wtime

      timetemp = wtime()

      if (np==0) return

      wpp = 1.

      oneondt = 1./dvnz(dt)

c       --- loop over particles
      do ip=1,np
        z_local = zp(ip)!-zbeam
        izl  = 1+int((z_local - tloc_zmin)*tloc_dzi)

c       --- cycle if particle out of zone of calculation
        if(izl<1 .or. izl>nztlocator) cycle

c       --- loop over temperature slices
        do its = 1, ntl(izl)
          izt = tslice_locator(izl,its)

          if(l_temp_collapseinz) then
c           --- collapse slice in z, i.e. align Z-locations of particles using
c           --- current velocity and from previous time step
            zc = 0.5*(tslicezmin(izt)+tslicezmax(izt))
            vzi = 1./(uzp(ip)*gaminv(ip)+SMALLPOS)
            dti  = (zbeam+zc-zp(ip))*vzi
            xpt  = xp(ip) + uxp(ip)*dti*gaminv(ip)
            ypt  = yp(ip) + uyp(ip)*dti*gaminv(ip)
            zpt  = z_local + uzp(ip)*dti*gaminv(ip)
            vxpt = (uxp(ip)*(1. + dti*oneondt) - uxpo(ip)*dti*oneondt) * gaminv(ip)
            vypt = (uyp(ip)*(1. + dti*oneondt) - uypo(ip)*dti*oneondt) * gaminv(ip)
            vzpt = (uzp(ip)*(1. + dti*oneondt) - uzpo(ip)*dti*oneondt) * gaminv(ip)
          else
            xpt  = xp(ip)
            ypt  = yp(ip)
            zpt  = z_local
            vxpt = uxp(ip) * gaminv(ip)
            vypt = uyp(ip) * gaminv(ip)
            vzpt = uzp(ip) * gaminv(ip)
          end if
          if(lrtheta) then
            rpt = sqrt(xpt**2+ypt**2)
            if(rpt>SMALLPOS) then
              cost = xpt/rpt
              sint = ypt/rpt
              xpt = vxpt
              vxpt = vxpt*cost + vypt*sint
              vypt = -xpt*sint + vypt*cost
            else
              vxpt = sqrt(vxpt**2+vypt**2)
              vypt = 0.
            end if
            xpt = rpt
            ypt = 0.
          else
            if(l2symtry) ypt = abs(ypt)
            if(l4symtry) then
              xpt=abs(xpt)
              ypt=abs(ypt)
            end if
          end if

c         --- cycle if particle not in slice
          if(lrtheta) then
            if(xpt>tslicexmax(izt) .or.
     &         zpt<=tslicezmin(izt) .or. zpt>tslicezmax(izt)) cycle
          else
            if(xpt<=tslicexmin(izt) .or. xpt>tslicexmax(izt) .or.
     &         ypt<=tsliceymin(izt) .or. ypt>tsliceymax(izt) .or.
     &         zpt<=tslicezmin(izt) .or. zpt>tslicezmax(izt)) cycle
          end if

c         --- compute location in arrays
          xt = (xpt-tslicexmin(izt))*dxti(izt)
          yt = (ypt-tsliceymin(izt))*dyti(izt)
          ixt = max(0,min(nxtslices-1,int(xt)))
          iyt = max(0,min(nytslices-1,int(yt)))

          if(lw) wpp = wp(ip)
          wt=wpp

c         --- deposit data in arrays
          pnumt  (ixt,iyt,izt) = pnumt  (ixt,iyt,izt) + 1.
          pnumtw (ixt,iyt,izt) = pnumtw (ixt,iyt,izt) + wt

          vxbart  (ixt,iyt,izt) = vxbart  (ixt,iyt,izt) + wt * vxpt
          vybart  (ixt,iyt,izt) = vybart  (ixt,iyt,izt) + wt * vypt
          vzbart  (ixt,iyt,izt) = vzbart  (ixt,iyt,izt) + wt * vzpt

          vxsqbart(ixt,iyt,izt) = vxsqbart(ixt,iyt,izt) + wt * vxpt**2
          vysqbart(ixt,iyt,izt) = vysqbart(ixt,iyt,izt) + wt * vypt**2
          vzsqbart(ixt,iyt,izt) = vzsqbart(ixt,iyt,izt) + wt * vzpt**2

          if (lrelativ) then
            ke = (1./gaminv(ip)-1.)
          else
            ke = vxpt*vxpt + vypt*vypt + vzpt*vzpt
          end if
          kebart  (ixt,iyt,izt) = kebart  (ixt,iyt,izt) + wt * ke
          kesqbart(ixt,iyt,izt) = kesqbart(ixt,iyt,izt) + wt * ke**2

          if(l_temp_rmcorrelations .or. l_temp_rmcrosscorrelations) then
            xbart   (ixt,iyt,izt) = xbart   (ixt,iyt,izt) + wt * xpt
            ybart   (ixt,iyt,izt) = ybart   (ixt,iyt,izt) + wt * ypt
            zbart   (ixt,iyt,izt) = zbart   (ixt,iyt,izt) + wt * zpt

            xsqbart (ixt,iyt,izt) = xsqbart (ixt,iyt,izt) + wt * xpt**2
            ysqbart (ixt,iyt,izt) = ysqbart (ixt,iyt,izt) + wt * ypt**2
            zsqbart (ixt,iyt,izt) = zsqbart (ixt,iyt,izt) + wt * zpt**2
          endif

          if(l_temp_rmcorrelations .or. l_temp_rmcrosscorrelations) then
            xvxbart (ixt,iyt,izt) = xvxbart (ixt,iyt,izt) + wt * xpt * vxpt
            yvybart (ixt,iyt,izt) = yvybart (ixt,iyt,izt) + wt * ypt * vypt
            zvzbart (ixt,iyt,izt) = zvzbart (ixt,iyt,izt) + wt * zpt * vzpt
            xkebart (ixt,iyt,izt) = xkebart (ixt,iyt,izt) + wt * xpt * ke
            ykebart (ixt,iyt,izt) = ykebart (ixt,iyt,izt) + wt * ypt * ke
            zkebart (ixt,iyt,izt) = zkebart (ixt,iyt,izt) + wt * zpt * ke
          endif

          if(l_temp_rmcrosscorrelations) then
            xybart (ixt,iyt,izt) = xybart (ixt,iyt,izt) + wt * xpt * ypt
            xzbart (ixt,iyt,izt) = xzbart (ixt,iyt,izt) + wt * xpt * zpt
            yzbart (ixt,iyt,izt) = yzbart (ixt,iyt,izt) + wt * ypt * zpt
            xvybart (ixt,iyt,izt) = xvybart (ixt,iyt,izt) + wt * xpt * vypt
            xvzbart (ixt,iyt,izt) = xvzbart (ixt,iyt,izt) + wt * xpt * vzpt
            yvxbart (ixt,iyt,izt) = yvxbart (ixt,iyt,izt) + wt * ypt * vxpt
            yvzbart (ixt,iyt,izt) = yvzbart (ixt,iyt,izt) + wt * ypt * vzpt
            zvxbart (ixt,iyt,izt) = zvxbart (ixt,iyt,izt) + wt * zpt * vxpt
            zvybart (ixt,iyt,izt) = zvybart (ixt,iyt,izt) + wt * zpt * vypt
          endif

        enddo
      enddo

      temperaturestime = temperaturestime + (wtime() - timetemp)
      return
      end
c===========================================================================
      subroutine finalize_temperature(is,m,lrtheta)
      use Constant
      use Temperatures
      use Timers
      use Beam_acc,only:lrelativ
      integer(ISZ):: is, ixt, iyt, izt
      real(kind=8):: m,timetemp,wtime
      real(kind=8):: pnumi,tfact,tottmp,kevar
      real(kind=8):: xvar,vxvar,xvxbar,xvybar,xvzbar,xybar,xkebar
      real(kind=8):: yvar,vyvar,yvxbar,yvybar,yvzbar,xzbar,ykebar
      real(kind=8):: zvar,vzvar,zvxbar,zvybar,zvzbar,yzbar,zkebar
      real(kind=8):: a1,a2,a3,b1,b2,b3,c1,c2,c3,d1x,d2x,d3x,d1y,d2y,d3y,d1z,d2z,d3z,d1ke,d2ke,d3ke,d
      real(kind=8)::xvxslope,xvyslope,xvzslope,xkeslope
      real(kind=8)::yvxslope,yvyslope,yvzslope,ykeslope
      real(kind=8)::zvxslope,zvyslope,zvzslope,zkeslope
      logical(ISZ)::lrtheta
      timetemp = wtime()

#ifdef MPIPARALLEL
c     --- For slave, call routine which sums moments over processors.
      call parallel_sum_temperature
#endif

c     --- set multiplying factor for proper units
      if(t_units==evolt)  tfact = 0.5*m/echarge
      if(t_units==joule)  tfact = 0.5*m
      if(t_units==kelvin) tfact = 0.5*m/boltzmann

c     --- Complete the calculation of temperatures: divide by particle number
      do izt = 1, nztslices
        do iyt = 0, nytslices
          do ixt = 0, nxtslices
            pnumi = 1./(pnumtw(ixt,iyt,izt)+SMALLPOS)

c           --- Compute averages
            vxbart  (ixt,iyt,izt) = vxbart  (ixt,iyt,izt) * pnumi
            vybart  (ixt,iyt,izt) = vybart  (ixt,iyt,izt) * pnumi
            vzbart  (ixt,iyt,izt) = vzbart  (ixt,iyt,izt) * pnumi

            vxsqbart(ixt,iyt,izt) = vxsqbart(ixt,iyt,izt) * pnumi
            vysqbart(ixt,iyt,izt) = vysqbart(ixt,iyt,izt) * pnumi
            vzsqbart(ixt,iyt,izt) = vzsqbart(ixt,iyt,izt) * pnumi

            kebart  (ixt,iyt,izt) = kebart  (ixt,iyt,izt) * pnumi
            kesqbart(ixt,iyt,izt) = kesqbart(ixt,iyt,izt) * pnumi

c           --- Compute second order moments with averages subtracted
            vxvar = vxsqbart(ixt,iyt,izt) - vxbart(ixt,iyt,izt)**2
            vyvar = vysqbart(ixt,iyt,izt) - vybart(ixt,iyt,izt)**2
            vzvar = vzsqbart(ixt,iyt,izt) - vzbart(ixt,iyt,izt)**2
            kevar = kesqbart(ixt,iyt,izt) - kebart(ixt,iyt,izt)**2

            if (pnumt(ixt,iyt,izt)<=5.) then
              tempx(ixt,iyt,izt,is) = 0.
              tempy(ixt,iyt,izt,is) = 0.
              tempz(ixt,iyt,izt,is) = 0.
              dke(ixt,iyt,izt,is) = 0.
              cycle
            end if

c           --- Compute linear correlations
            if(l_temp_rmcorrelations .or. l_temp_rmcorrelations) then
              xbart   (ixt,iyt,izt) = xbart   (ixt,iyt,izt) * pnumi
              ybart   (ixt,iyt,izt) = ybart   (ixt,iyt,izt) * pnumi
              zbart   (ixt,iyt,izt) = zbart   (ixt,iyt,izt) * pnumi

              xsqbart (ixt,iyt,izt) = xsqbart (ixt,iyt,izt) * pnumi
              ysqbart (ixt,iyt,izt) = ysqbart (ixt,iyt,izt) * pnumi
              zsqbart (ixt,iyt,izt) = zsqbart (ixt,iyt,izt) * pnumi

              xvar  = xsqbart (ixt,iyt,izt) - xbart (ixt,iyt,izt)**2
              yvar  = ysqbart (ixt,iyt,izt) - ybart (ixt,iyt,izt)**2
              zvar  = zsqbart (ixt,iyt,izt) - zbart (ixt,iyt,izt)**2

              xvxslope = 0.
              yvxslope = 0.
              zvxslope = 0.
              xvyslope = 0.
              yvyslope = 0.
              zvyslope = 0.
              xvzslope = 0.
              yvzslope = 0.
              zvzslope = 0.
              xkeslope = 0.
              ykeslope = 0.
              zkeslope = 0.

              if (.not. l_temp_rmcrosscorrelations) then
                xvxbart (ixt,iyt,izt) = xvxbart (ixt,iyt,izt) * pnumi
                yvybart (ixt,iyt,izt) = yvybart (ixt,iyt,izt) * pnumi
                zvzbart (ixt,iyt,izt) = zvzbart (ixt,iyt,izt) * pnumi
                xvxbar  = xvxbart (ixt,iyt,izt) - xbart (ixt,iyt,izt)*vxbart(ixt,iyt,izt)
                yvybar  = yvybart (ixt,iyt,izt) - ybart (ixt,iyt,izt)*vybart(ixt,iyt,izt)
                zvzbar  = zvzbart (ixt,iyt,izt) - zbart (ixt,iyt,izt)*vzbart(ixt,iyt,izt)
                if(abs(xvar)>SMALLPOS) vxvar = vxvar - xvxbar**2/xvar
                if(abs(yvar)>SMALLPOS) vyvar = vyvar - yvybar**2/yvar
                if(abs(zvar)>SMALLPOS) vzvar = vzvar - zvzbar**2/zvar
              else ! l_temp_rmcrosscorrelations is true
                xvxbart (ixt,iyt,izt) = xvxbart (ixt,iyt,izt) * pnumi
                xvybart (ixt,iyt,izt) = xvybart (ixt,iyt,izt) * pnumi
                xvzbart (ixt,iyt,izt) = xvzbart (ixt,iyt,izt) * pnumi
                yvxbart (ixt,iyt,izt) = yvxbart (ixt,iyt,izt) * pnumi
                yvybart (ixt,iyt,izt) = yvybart (ixt,iyt,izt) * pnumi
                yvzbart (ixt,iyt,izt) = yvzbart (ixt,iyt,izt) * pnumi
                zvxbart (ixt,iyt,izt) = zvxbart (ixt,iyt,izt) * pnumi
                zvybart (ixt,iyt,izt) = zvybart (ixt,iyt,izt) * pnumi
                zvzbart (ixt,iyt,izt) = zvzbart (ixt,iyt,izt) * pnumi
                xkebart (ixt,iyt,izt) = xkebart (ixt,iyt,izt) * pnumi
                ykebart (ixt,iyt,izt) = ykebart (ixt,iyt,izt) * pnumi
                zkebart (ixt,iyt,izt) = zkebart (ixt,iyt,izt) * pnumi
                xybart (ixt,iyt,izt) = xybart (ixt,iyt,izt) * pnumi
                xzbart (ixt,iyt,izt) = xzbart (ixt,iyt,izt) * pnumi
                yzbart (ixt,iyt,izt) = yzbart (ixt,iyt,izt) * pnumi
                xvxbar  = xvxbart (ixt,iyt,izt) - xbart (ixt,iyt,izt)*vxbart(ixt,iyt,izt)
                xvybar  = xvybart (ixt,iyt,izt) - xbart (ixt,iyt,izt)*vybart(ixt,iyt,izt)
                xvzbar  = xvzbart (ixt,iyt,izt) - xbart (ixt,iyt,izt)*vzbart(ixt,iyt,izt)
                yvxbar  = yvxbart (ixt,iyt,izt) - ybart (ixt,iyt,izt)*vxbart(ixt,iyt,izt)
                yvybar  = yvybart (ixt,iyt,izt) - ybart (ixt,iyt,izt)*vybart(ixt,iyt,izt)
                yvzbar  = yvzbart (ixt,iyt,izt) - ybart (ixt,iyt,izt)*vzbart(ixt,iyt,izt)
                zvxbar  = zvxbart (ixt,iyt,izt) - zbart (ixt,iyt,izt)*vxbart(ixt,iyt,izt)
                zvybar  = zvybart (ixt,iyt,izt) - zbart (ixt,iyt,izt)*vybart(ixt,iyt,izt)
                zvzbar  = zvzbart (ixt,iyt,izt) - zbart (ixt,iyt,izt)*vzbart(ixt,iyt,izt)
                xkebar  = xkebart (ixt,iyt,izt) - xbart (ixt,iyt,izt)*kebart(ixt,iyt,izt)
                ykebar  = ykebart (ixt,iyt,izt) - ybart (ixt,iyt,izt)*kebart(ixt,iyt,izt)
                zkebar  = zkebart (ixt,iyt,izt) - zbart (ixt,iyt,izt)*kebart(ixt,iyt,izt)
                xybar  = xybart (ixt,iyt,izt) - xbart (ixt,iyt,izt)*ybart(ixt,iyt,izt)
                xzbar  = xzbart (ixt,iyt,izt) - xbart (ixt,iyt,izt)*zbart(ixt,iyt,izt)
                yzbar  = yzbart (ixt,iyt,izt) - ybart (ixt,iyt,izt)*zbart(ixt,iyt,izt)
                ! --- sets coefficients of linear system
                a1 = xvar;  b1 = xybar; c1 = xzbar; d1x = xvxbar; d1y = xvybar; d1z = xvzbar; d1ke = xkebar
                a2 = xybar; b2 = yvar;  c2 = yzbar; d2x = yvxbar; d2y = yvybar; d2z = yvzbar; d2ke = ykebar
                a3 = xzbar; b3 = yzbar; c3 = zvar;  d3x = zvxbar; d3y = zvybar; d3z = zvzbar; d3ke = zkebar
                if (lrtheta) then
                 d = (a1*c3-a3*c1)
                 if (abs(d)>SMALLPOS) then
                  d=1./d
                  xvxslope = (c1*d3x-d1x*c3)*d
                  yvxslope = 0.
                  zvxslope = -(a3*xvxslope+d3x)/c3
                  xvyslope = (c1*d3y-d1y*c3)*d
                  yvyslope = 0.
                  zvyslope = -(a3*xvyslope+d3y)/c3
                  xvzslope = (c1*d3z-d1z*c3)*d
                  yvzslope = 0.
                  zvzslope = -(a3*xvzslope+d3z)/c3
                  xkeslope = (c1*d3ke-d1ke*c3)*d
                  ykeslope = 0.
                  zkeslope = -(a3*xkeslope+d3ke)/c3
                 end if
                else
                 d = ((a1*c2-a2*c1)*(b2*c3-b3*c2)-(a2*c3-a3*c2)*(b1*c2-b2*c1))
                 if (abs(d)>SMALLPOS) then
                  d=1./d
                  xvxslope = ((d2x*c3-d3x*c2)*(b1*c2-b2*c1)-(d1x*c2-d2x*c1)*(b2*c3-b3*c2))*d
                  if (abs(b3*c2-b2*c3)>SMALLPOS) yvxslope = ((a2*c3-a3*c2)*xvxslope+d2x*c3-d3x*c2)/(b3*c2-b2*c3)
                  if (abs(c1)>SMALLPOS) zvxslope = -(a1*xvxslope+b1*yvxslope+d1x)/c1
                  xvyslope = ((d2y*c3-d3y*c2)*(b1*c2-b2*c1)-(d1y*c2-d2y*c1)*(b2*c3-b3*c2))*d
                  if (abs(b3*c2-b2*c3)>SMALLPOS) yvyslope = ((a2*c3-a3*c2)*xvyslope+d2y*c3-d3y*c2)/(b3*c2-b2*c3)
                  if (abs(c1)>SMALLPOS) zvyslope = -(a1*xvyslope+b1*yvyslope+d1y)/c1
                  xvzslope = ((d2z*c3-d3z*c2)*(b1*c2-b2*c1)-(d1z*c2-d2z*c1)*(b2*c3-b3*c2))*d
                  if (abs(b3*c2-b2*c3)>SMALLPOS) yvzslope = ((a2*c3-a3*c2)*xvzslope+d2z*c3-d3z*c2)/(b3*c2-b2*c3)
                  if (abs(c1)>SMALLPOS) zvzslope = -(a1*xvzslope+b1*yvzslope+d1z)/c1
                  xkeslope = ((d2ke*c3-d3ke*c2)*(b1*c2-b2*c1)-(d1ke*c2-d2ke*c1)*(b2*c3-b3*c2))*d
                  if (abs(b3*c2-b2*c3)>SMALLPOS) ykeslope = ((a2*c3-a3*c2)*xkeslope+d2ke*c3-d3ke*c2)/(b3*c2-b2*c3)
                  if (abs(c1)>SMALLPOS) zkeslope = -(a1*xkeslope+b1*ykeslope+d1ke)/c1
                 end if
                end if
                ! --- remove correlations on vx, vy and vz
                vxvar = vxvar - xvxslope**2*xvar - yvxslope**2*yvar - zvxslope**2*zvar
     &                        - 2*xvzslope*yvxslope*xybar
     &                        - 2*xvzslope*zvxslope*xzbar
     &                        - 2*yvzslope*zvxslope*yzbar
                vyvar = vyvar - xvyslope**2*xvar - yvyslope**2*yvar - zvyslope**2*zvar
     &                        - 2*xvzslope*yvyslope*xybar
     &                        - 2*xvzslope*zvyslope*xzbar
     &                        - 2*yvzslope*zvyslope*yzbar
                vzvar = vzvar - xvzslope**2*xvar - yvzslope**2*yvar - zvzslope**2*zvar
     &                        - 2*xvzslope*yvzslope*xybar
     &                        - 2*xvzslope*zvzslope*xzbar
     &                        - 2*yvzslope*zvzslope*yzbar
                vxvar = abs(vxvar)
                vyvar = abs(vyvar)
                vzvar = abs(vzvar)
                kevar = kevar - xkeslope**2*xvar - ykeslope**2*yvar - zkeslope**2*zvar
     &                        - 2*xkeslope*ykeslope*xybar
     &                        - 2*xkeslope*zkeslope*xzbar
     &                        - 2*ykeslope*zkeslope*yzbar
              end if
            end if

c           --- Compute temperatures
            tempx(ixt,iyt,izt,is) = tfact*vxvar
            tempy(ixt,iyt,izt,is) = tfact*vyvar
            tempz(ixt,iyt,izt,is) = tfact*vzvar

c           --- Compute energy spread
            dke(ixt,iyt,izt,is)   = tfact*sqrt(kevar)
            if (lrelativ) dke(ixt,iyt,izt,is) = dke(ixt,iyt,izt,is) * 2.*clight**2

          enddo
        enddo
        tottmp = sum(pnumtw(:,:,izt))
        if(tottmp>SMALLPOS) then
          tempxz(izt,is) = sum(tempx(:,:,izt,is)*pnumtw(:,:,izt))/tottmp
          tempyz(izt,is) = sum(tempy(:,:,izt,is)*pnumtw(:,:,izt))/tottmp
          tempzz(izt,is) = sum(tempz(:,:,izt,is)*pnumtw(:,:,izt))/tottmp
        end if
      enddo
      temperaturestime = temperaturestime + (wtime() - timetemp)
      return
      end
c=============================================================================
      subroutine gett(is,lrtheta,l2symtry,l4symtry)
c     Compute temperature in Z-slices for species 'is' on a 3-D grid:
c       - the slices can have any position and thickness but cannot overlap,
c       - the min and max of each slice in x, y and z are given respectively in the arrays
c         tslicexmin, tslicexmax, tsliceymin, tsliceymax, tslicezmin and tslicezmax
c         (the reason for having slices with different thickness and dimensions is
c          to allow the temperature measurement to the shape of a distribution, like
c          for example a beam extending over several quadrupoles and accelerating gaps),
c       - the x, y and z temperatures are given in the arrays tempx, tempy and tempz,
c         while averages in each slice are given in the arrays tempxz, tempyz, tempzz.
c       - the calculation is done in three parts:
c         o reset_temperature: zero out all moments
c         o accumulate_temperature: accumulate moments from particles
c         o finalize_temperature: compute final quantities
c       - l_temp_collapseinz=.true.: collapse slices in z, i.e. align Z-locations
c                                    of particles using current velocity (uxp,uyp,yzp) and
c                                    velocity from previous time step (uxpo,uypo,uzpo)
c       - lrtheta=.true.: radial and azimuthal are computed in place of X and Y.
c       - the default temperature unit is in electron-volt. To select the units,
c         set the variable t_unit to the default integers evolt, joule or kelvin.
        use InGen
        use Particles,Only: pgroup,wpid
        integer(ISZ):: is, ipmin, itask, i1, i2
        logical(ISZ):: lrtheta,l2symtry,l4symtry

        ipmin = pgroup%ins(is)
        call reset_temperature(is)
        if (pgroup%nps(is)>0) then
          i1 = ipmin
          i2 = ipmin + pgroup%nps(is) - 1
          if(wpid>0) then
            call accumulate_temperature(pgroup%nps(is),
     &                pgroup%xp(i1:i2),pgroup%yp(i1:i2),pgroup%zp(i1:i2),
     &                pgroup%uxp(i1:i2),pgroup%uyp(i1:i2),pgroup%uzp(i1:i2),
     &                pgroup%gaminv(i1:i2),pgroup%sw(is),dt,
     &                pgroup%uxp(i1:i2),pgroup%uyp(i1:i2),pgroup%uzp(i1:i2),is,
     &                pgroup%pid(i1:i2,wpid),.true.,lrtheta,l2symtry,l4symtry)
          else
            call accumulate_temperature(pgroup%nps(is),
     &                pgroup%xp(i1:i2),pgroup%yp(i1:i2),pgroup%zp(i1:i2),
     &                pgroup%uxp(i1:i2),pgroup%uyp(i1:i2),pgroup%uzp(i1:i2),
     &                pgroup%gaminv(i1:i2),pgroup%sw(is),dt,
     &                pgroup%uxp(i1:i2),pgroup%uyp(i1:i2),pgroup%uzp(i1:i2),is,
     &                pgroup%xp(i1:i2),.false.,lrtheta,l2symtry,l4symtry)
          end if
        end if
        call finalize_temperature(is,pgroup%sm(is),lrtheta)
        return
      end
c=============================================================================
      subroutine setregulartgrid(nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,dz,nzloc,lcollapse,lcorrel,lcrosscorrel)
c       Setup regular grid for temperature calculation. The temperature is calculated in nz slices
c       evenly spaced of thickness dz. The spacing between each slice is given by (zmax-zmin)/(nz-1).
c       In each slice, the temperature will be computed on a nx+1*ny+1 grid of size (xmin,xmax,ymin,ymax).
        use Temperatures
        use InPart

        integer(ISZ) :: nx, ny, nz, nzloc
        logical(ISZ) :: lcollapse, lcorrel, lcrosscorrel
        real(kind=8) :: xmin, xmax, ymin, ymax, zmin, zmax, dz

        integer(ISZ) :: i

        nstemp = ns

        nxtslices = nx
        nytslices = ny
        nztslices = nz
        if(lcorrel) then
          l_temp_rmcorrelations = .true.
          nxtslicesc = nx
          nytslicesc = ny
          nztslicesc = nz
        else
          l_temp_rmcorrelations = .false.
        end if
        if(lcrosscorrel) then
          l_temp_rmcrosscorrelations = .true.
          nxtslicescc = nx
          nytslicescc = ny
          nztslicescc = nz
        else
          l_temp_rmcrosscorrelations = .false.
        end if
        nztlocator = nzloc
        call gchange("Temperatures",0)
        tslicexmin = xmin
        tslicexmax = xmax
        tsliceymin = ymin
        tsliceymax = ymax
        if (xmax .ne. xmin) then
          dxti = nx/(xmax-xmin)
        else
          dxti = 0.
        endif
        if (ymax .ne. ymin) then
          dyti = ny/(ymax-ymin)
        else
          dyti = 0.
        endif
        if(nztslices>1) then
          do i = 1, nztslices
            tslicezmin(i) = zmin+(i-1)*(zmax-zmin)/(nztslices-1)-0.5*dz
            tslicezmax(i) = tslicezmin(i) + dz
          end do
        else
          tslicezmin = zmin-0.5*dz
          tslicezmax = zmin+0.5*dz
        end if
        tloc_dzi = nzloc/(tslicezmax(nztslices)-tslicezmin(1))
        tloc_zmin = tslicezmin(1)
        tloc_zmax = tslicezmax(nztslices)
        call set_tslice_locator()

        return
      end
c=============================================================================
      subroutine set_tslice_locator()
        use Temperatures
        integer(ISZ) :: i, ii, izmin, izmax

        integer(ISZ), allocatable :: ntsloc(:,:)
        integer(ISZ):: allocerror

        allocate(ntsloc(nztlocator,nztslices),stat=allocerror)
        if (allocerror /= 0) then
          print*,"set_tslice_locator: allocation error ",allocerror,
     &           ": could not allocate ntsloc to shape ",nztlocator,nztslices
          call kaboom("set_tslice_locator: allocation error")
          return
        endif

        ntsloc = 0
        ntl = 0
        do i = 1, nztslices
          izmin = min(nztlocator,max(1,1+int((tslicezmin(i)-tloc_zmin)*tloc_dzi)))
          izmax = min(nztlocator,max(1,1+int((tslicezmax(i)-tloc_zmin)*tloc_dzi)))
          do ii = izmin,izmax
            ntl(ii) = ntl(ii)+1
            ntsloc(ii,ntl(ii)) = i
          end do
        end do
        ntlmax = maxval(ntl)
        call gchange("Temperatures",0)
        tslice_locator = ntsloc(:,:ntlmax)
c there is a problem here with the intel compiler
        deallocate(ntsloc)

        return
      end
c=============================================================================