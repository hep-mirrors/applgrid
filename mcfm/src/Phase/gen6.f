      subroutine gen6(r,q,wt6,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'phasemin.f'
      include 'process.f'
      integer nu
      double precision r(mxdim)
      double precision wt6,q(mxpart,4)
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4)
      double precision sqrts,y,pswt,xjac,xx(2),tau
      common/energy/sqrts
      common/x1x2/xx
      data p3/0d0,0d0,0d0,0d0/

      wt6=0d0

      tau=dexp(dlog(taumin)*r(9))
      y=0.5d0*dlog(tau)*(1d0-2d0*r(10))
      xjac=dlog(taumin)*tau*dlog(tau)

      xx(1)=dsqrt(tau)*dexp(+y)
      xx(2)=dsqrt(tau)*dexp(-y)

c--- phase space volume only checked for x1=x2=1
      if ((case .eq. 'vlchwg') .or. (case .eq. 'vlchwh')) then
        xx(1)=1d0
        xx(2)=1d0
        xjac=1d0
      endif

c---if x's out of normal range alternative return
      if   ((xx(1) .gt. 1d0)
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) return 1

      p1(4)=-xx(1)*sqrts*half
      p1(1)=zip
      p1(2)=zip
      p1(3)=-xx(1)*sqrts*half

      p2(4)=-xx(2)*sqrts*half
      p2(1)=zip
      p2(2)=zip
      p2(3)=+xx(2)*sqrts*half

      if     ((case .eq. 'W_twdk') .or. (case .eq. 'Wtbwdk')
     .   .or. (case .eq. 'W_cwdk') .or. (case .eq. 'vlchwh')) then
c--- W+t process, radiation in production
        call phase6a(r,p1,p2,p3,p4,p5,p6,p7,p8,pswt,*999)
      elseif ((case .eq. 'Wtdkay') .or. (case .eq. 'vlchwg')) then
c--- W+t process, radiation in decay
        call phase6b(r,p1,p2,p3,p4,p5,p6,p7,p8,pswt,*999)
      else
c--- generic case 
        call phase6(r,p1,p2,p3,p4,p5,p6,p7,p8,pswt,*999)
      endif
      
      do nu=1,4
        q(1,nu)=p1(nu)
        q(2,nu)=p2(nu)
        q(3,nu)=p3(nu)
        q(4,nu)=p4(nu)
        q(5,nu)=p5(nu)
        q(6,nu)=p6(nu)
        q(7,nu)=p7(nu)
        q(8,nu)=p8(nu)
      enddo 
      
      wt6=xjac*pswt

      if (debug) write(6,*) 'wt6 in gen6',wt6

      return

 999  return 1
      end

