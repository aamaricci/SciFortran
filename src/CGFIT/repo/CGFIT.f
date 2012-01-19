!#####################################################################
!     PROGRAM  : FITGREEN
!     TYPE     : subroutine
!     PURPOSE  : fit given Green function with 
!     \Delta=Sum_i=1^nsite b(i)^2/(z-a(i))
!     using conjugate gradient
!     COMMENTS : 
!#####################################################################
!     ----------------------------------------------------------------
!     This subroutine fits a given GF g at frequencies z=i  w_n
!     to obtain a GF of the form
!     G(z) =Sum_i=1^nsite b(i)^2/(z-a(i))
!     using the conjugate gradient method. The number of iterations
!     used is given by iter, the chi^2 of the fit is given by chisq.
!     The number of points fitted is given by n.
!     Notice that the initial guesses for a's and b's should be
!     different (scattered) even in the case of half filling.
!     At half filling choosing them symmetrically is useful due
!     to thes sensitivity of the CG method on the initial guess.
!     THe fitted coefficients are in the vectors a(nsite) and
!     b(nsite).
!     na is the number of adjustable parameters.
      
      subroutine fitgreen(nf,z,g,nsite,agues,bgues,aa,bb,aamin)

      integer nf,  iter, nsite, i, n, na, nff,NL
      integer LL
      parameter(LL=2**14)       !change here AND in the other 2places
      parameter(nmax = 100)
      double complex g(*), z(*)
      double precision aa(nsite),bb(nsite),agues(nsite),bgues(nsite)
      double precision xome(LL), gre(LL), gim(LL), ftol
      double precision a(nmax), chi, aamin, min
      
      common/vars/xome,gre,gim,nff,na
      common/vars2/min

      min=aamin/2

      na=2*nsite
      nff = nf                  !4
      if(nff.gt.1024)nff=1024   !fix an upper bound
      print*, 'n_fit=', nff

      do i=1,nf
         xome(i)= dimag(z(i))
         gre(i) =  dble(g(i))
         gim(i) = dimag(g(i))
      enddo

      do i=1,nsite
         a(i)       = agues(i)
         a(i+nsite) = bgues(i)
      enddo

      ftol = 1.e-12             !.00000000001d0
      call conjugate(a,na,ftol,iter,chi)

      do i=1,nsite
         aa(i) = a(i)
         bb(i) = a(i+nsite)
      enddo
      print*,'chi^2,iter=',chi,iter
      print*,'------------------------------------'
      return
      end
!     ================================================================


!#####################################################################
!     PROGRAM  : CONJUGATE
!     TYPE     : subroutine
!     PURPOSE  : Minimize the Chi^2 distance using conjugate gradient
!     COMMENTS : 
!#####################################################################
      subroutine conjugate(p,na,ftol,iter,fret)
!     Adapted by FRPRM subroutine from NumRec (10.6)
C     Given a starting point P that is a vector of length N, 
!     the Fletcher-Reeves-Polak-Ribiere minimisation is performed 
!     n a functin FUNC,using its gradient as calculated by a 
!     routine DFUNC. The convergence tolerance on the function 
!     value is input as FTOL.  
!     Returned quantities are: 
!     - P (the location of the minimum), 
!     - ITER (the number of iterations that were performed), 
!     - FRET (the minimum value of the function). 
!     The routine LINMIN is called to perform line minimisations.
     
!     Minimisation routines: DFPMIN, LINMIN, MNBRAK, BRENT and F1DIM
!     come from Numerical Recipes.

      integer itmax, nmax, iter
      double precision fret, ftol, p(na), eps, func
      external dfunc, func, linmin
      parameter(nmax = 100, itmax = 10000, eps = 1.d-6)
      integer its, j
      double precision dgg, fp, gam, gg, g(nmax), h(nmax), xi(nmax)

      integer na

      fp = func(p)
      call dfunc(p, xi)
      do i=1,na
         g(i) = -xi(i)
         h(i) = g(i)
         xi(i)= h(i)
      enddo

      do its=1, itmax           !Loop over the iterations
         iter = its
         call linmin(p,xi,na,fret)
!     return if converged to a minimum
         if (2.*dabs(fret-fp).le.ftol*(dabs(fret)+dabs(fp)+eps)) return 
!     fp=fret
         fp = func(p)
         call dfunc(p, xi)
         
         gg = 0.d0
         dgg= 0.d0
         do j=1,na
            gg = gg+g(j)**2
!     dgg=dgg+xi(j)**2    !Fletcher-Reeves.
            dgg=dgg+(xi(j)+g(j))*xi(j) !Polak-Ribiere
         enddo
!     Unlikely if the gradient is 0.d0 then we're already done.
         if (gg.eq.0.d0) return 
         gam = dgg/gg
         do j=1,na
            g(j) =-xi(j)
            h(j) =g(j)+gam*h(j)
            xi(j)=h(j)
         enddo
      enddo
      print*, 'Too many iteration in CG'
      return
      end
!======================================================================


!########################################################################
!     PROGRAM  : FUNC
!     TYPE     : subroutine
!     PURPOSE  : Obtain Chi^2 distance
!     COMMENTS : 
!########################################################################
      double precision function func(a)
      integer  nmax, i, n, na, nff
      parameter(nmax = 100)
      integer LL
      parameter(LL=2**14)
      double precision a(nmax), yrep, yimp, min,om
      double precision gre(LL), gim(LL), xome(LL)
      complex*16 yand,gand
      common/vars/ xome, gre, gim, nff, na
      common/vars2/ min

      func = 0.d0
      do i=1,nff
         om=1.d0                !(xome(i))!tanh(xome(i))
         if (dabs(xome(i)).gt.dabs(min)) then
            call getpacf(xome(i), a, na, yrep, yimp)
            func = func+((gre(i)-yrep)**2+(gim(i)-yimp)**2)/om
         endif
      enddo
      return
      end
!==================================================================


!########################################################################
!     PROGRAM  : DFUNC
!     TYPE     : subroutine
!     PURPOSE  : obtain the gradient of Chi^2 function
!     COMMENTS : 
!########################################################################
      subroutine dfunc(a, df)
      integer nmax, i, j, na, nff!, na
      parameter(nmax = 100)
      integer LL
      parameter(LL=2**14)
      double precision a(nmax), df(nmax), dyrep(nmax), dyimp(nmax),min
      double precision gre(LL), gim(LL), xome(LL), yrep, yimp,sgn,om
      real*8 reg,img,rsgn,isgn
      common/vars/ xome,gre,gim,nff,na
      common/vars2/ min

      do i=1,na
         df(i) = 0.d0
      enddo
      do i=1,nff
         om=1.d0                !(xome(i))!tanh(xome(i))
         if (dabs(xome(i)).gt.dabs(min)) then
            call  getpacf(xome(i), a, na, yrep, yimp)
            call getdpacf(xome(i), a, na, dyrep, dyimp)
            do j=1,na
               df(j)=df(j)+((gre(i)-yrep)*dyrep(j)+
     *              (gim(i)-yimp)*dyimp(j))/om
            enddo
         endif
      enddo
      do j=1,na
         df(j)=-2.d0*df(j)
      enddo
      return
      end
!=======================================================================


!########################################################################
!     PROGRAM  : GETPACF
!     TYPE     : subroutine
!     PURPOSE  : Construct the Anderson function to fit against given GF
!     COMMENTS : 
!########################################################################
      subroutine getpacf(x,a,na,yrep,yimp)
      integer na, i, nmax
      parameter(nmax=100)
      double precision x, yrep, yimp, a(nmax)
      double complex z, gz

      ns=na/2
      z = dcmplx(0.d0,x)
      gz =0.d0 
      do i=1,ns                 !na/2
         gz = gz + a(i+ns)**2/(z-a(i))
      enddo
      yrep = dble(gz)
      yimp = dimag(gz)
      return
      end
!====================================================================


!########################################################################
!     PROGRAM  : GETDPACF
!     TYPE     : subroutine
!     PURPOSE  : get the gradient of the Anderson gf
!     COMMENTS : 
!########################################################################
      subroutine getdpacf(x, a, na, dyrep, dyimp)
      integer na, i,  nmax, ns
      parameter(nmax=100)
      double precision x, a(nmax),  dyimp(nmax), dyrep(nmax)
      double complex z, dgz(nmax),gz

      z = dcmplx(0.,x)
      ns=na/2
      do i=1,ns
         dgz(i)    = a(i+ns)*a(i+ns)/(z-a(i))**2
         dgz(i+ns) =      2.d0*a(i+ns)/(z-a(i))
      enddo
      do i=1,na
         dyrep(i) =  dble(dgz(i))
         dyimp(i) = dimag(dgz(i))
      enddo
      return
      end
!======================================================================


!########################################################################
!     PROGRAM  : LINMIN
!     TYPE     : subroutine
!     PURPOSE  : Minimization routine
!     COMMENTS : 
!########################################################################
      subroutine linmin(p, xi, na, fret)
!     Given an N dimensional point P and an N dimensional direction 
!     XI, LINMIN
!     moves and resets P to where the function FUNC(P) takes on a minimum
!     along the direction XI from P, and replaces XI by the actual vector
!     displacement that P was moved.  Also returns FRET the value of 
!     FUNC at the returned location P.  
!     This is actually all accomplished by calling the routines 
!     MNBRAK and BRENT.
      integer na, nmax
      double precision fret, p(na), xi(na), tol
      parameter(nmax=100, tol=1.e-4)
      integer j, ncom
      double precision ax,bx,fa,fb,fx,xmin,xx,pcom(nmax),xicom(nmax)
      double precision dbrent
      common/f1com/ pcom, xicom, ncom
      external f1dim,df1dim, dbrent, mnbrak

!...set up the common block
      ncom = na
      do j=1,na
         pcom(j) = p(j)
         xicom(j)= xi(j)
      enddo
      ax = 0.
      xx = 1.
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=dbrent(ax,xx,bx,f1dim,df1dim,tol,xmin)
!...construct the vector results to return
      do j=1,na
         xi(j)=xmin*xi(j)
         p(j) =p(j)+xi(j)
      enddo
      return
      end
!=====================================================================


!########################################################################
!     PROGRAM  : F1DIM
!     TYPE     : subroutine
!     PURPOSE  : 
!     COMMENTS : 
!########################################################################
      double precision function f1dim(x)
      integer nmax
      double precision func, x
      parameter(nmax=100)
      integer j,ncom
      double precision pcom(nmax), xicom(nmax),xt(nmax)
      common/f1com/ pcom, xicom, ncom
      external func
      do j=1,ncom
         xt(j)=pcom(j)+x*xicom(j)
      enddo
      f1dim = func(xt)
      return
      end
!======================================================================


!########################################################################
!     PROGRAM  : DF1DIM
!     TYPE     : subroutine
!     PURPOSE  : 
!     COMMENTS : 
!########################################################################
      double precision function df1dim(x)
      integer nmax
      double precision x
      parameter(nmax=100)
      integer j,ncom
      double precision df(nmax),  pcom(nmax), xicom(nmax),xt(nmax)
      common/f1com/ pcom, xicom, ncom
      external dfunc
      do j=1,ncom
         xt(j)=pcom(j)+x*xicom(j)
      enddo
      call dfunc(xt, df)
      df1dim = 0.
      do j=1,ncom
         df1dim=df1dim+df(j)*xicom(j)
      enddo
      return
      end
!=====================================================================


!########################################################################
!     PROGRAM  : MNBRAK
!     TYPE     : subroutine
!     PURPOSE  : 
!     COMMENTS : 
!########################################################################
      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
!     Given a function FUNC, and given distinct initial points AX and BX,
!     this routine searches in the downhill direction (defined by the 
!     function as evaluated at the initial points) and returns new points
!     AX, BX, CX which bracket a minimum of the function.  
!     Also returned are the function values at the three points, 
!     FA, FB, and FC.
      double precision ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
!...The first parameter is the default ratio by which successive intervals
!   are magnified; the second is the maximum magnification allowed for a
!   parabolic-fit step
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.d-20)
      double precision  dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     if(fb.ge.fc)then
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(dabs(q-r),TINY),q-r))
         ulim=bx+GLIMIT*(cx-bx)
         if((bx-u)*(u-cx).gt.0.)then
            fu=func(u)
            if(fu.lt.fc)then
               ax=bx
               fa=fb
               bx=u
               fb=fu
               return
            else if(fu.gt.fb)then
               cx=u
               fc=fu
               return
            endif
            u=cx+GOLD*(cx-bx)
            fu=func(u)
         else if((cx-u)*(u-ulim).gt.0.)then
            fu=func(u)
            if(fu.lt.fc)then
               bx=cx
               cx=u
               u=cx+GOLD*(cx-bx)
               fb=fc
               fc=fu
               fu=func(u)
            endif
         else if((u-ulim)*(ulim-cx).ge.0.)then
            u=ulim
            fu=func(u)
         else
            u=cx+GOLD*(cx-bx)
            fu=func(u)
         endif
         ax=bx
         bx=cx
         cx=u
         fa=fb
         fb=fc
         fc=fu
         goto 1
      endif
      return
      END
!======================================================================


!########################################################################
!     PROGRAM  : 
!     TYPE     : subroutine
!     PURPOSE  : 
!     COMMENTS : 
!########################################################################
      double precision FUNCTION dbrent(ax,bx,cx,f,df,tol,xmin)
!     Given a function F, and given a bracketing triplet of 
!     abscissas AX, BX, CX (such that BX is between AX and CX, 
!     and F(BX) is less than both F(AX) and F(CX)), this routine 
!     isolates the minimum to a fractional precision of about TOL 
!     using Brent's method. The abscissa of the minimum is returned 
!     as XMIN, and the minimum function value is returned as BRENT, 
!     the returned function value.
      INTEGER ITMAX
      double precision ax,bx,cx,tol,xmin,df,f,ZEPS
      EXTERNAL df,f
      PARAMETER (ITMAX=1000,ZEPS=1.0d-10)
      INTEGER iter
      double precision a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw
      double precision fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm
      LOGICAL ok1,ok2
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f(x)
      fv=fx
      fw=fx
      dx=df(x)
      dv=dx
      dw=dx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*dabs(x)+ZEPS
        tol2=2.*tol1
        if(dabs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(dabs(e).gt.tol1) then
          d1=2.*(b-a)
          d2=d1
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0.).and.(dx*d1.le.0.)
          ok2=((a-u2)*(u2-b).gt.0.).and.(dx*d2.le.0.)
          olde=e
          e=d
          if(.not.(ok1.or.ok2))then
            goto 1
          else if (ok1.and.ok2)then
            if(dabs(d1).lt.dabs(d2))then
              d=d1
            else
              d=d2
            endif
          else if (ok1)then
            d=d1
          else
            d=d2
          endif
          if(dabs(d).gt.dabs(0.5*olde))goto 1
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(dx.ge.0.) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5*e
2       if(dabs(d).ge.tol1) then
          u=x+d
          fu=f(u)
        else
          u=x+sign(tol1,d)
          fu=f(u)
          if(fu.gt.fx)goto 3
        endif
        du=df(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
         if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
          endif
        endif
11    continue
      pause 'dbrent exceeded maximum iterations'
3     xmin=x
      dbrent=fx
      return
      END


