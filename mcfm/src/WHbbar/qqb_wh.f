      subroutine qqb_wh(p,msq)
c---Matrix element squared averaged over initial colors and spins
c for nwz=1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+b(p6)
c for nwz=-1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> e^-(p3)+nubar(p4)
c                           |
c                           ---> b(p5)+b(p6)
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ckm.f'
      integer j,k
      double precision p(mxpart,4)
      double precision s,prop,fac,qqbWH,qbqWH,s56
      double precision msq(-nf:nf,-nf:nf)

      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      s56=s(5,6)+2*mb**2
c---cut to ensure hard process
      if (
     .      (s(5,6) .lt. four*mbsq) 
     . .or. (s(1,5)*s(2,5)/s(1,2) .lt. mbsq) 
     . .or. (s(1,6)*s(2,6)/s(1,2) .lt. mbsq) ) return

c---calculate the 3 propagators
      prop=     ((s(1,2)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s56-hmass**2)**2+(hmass*hwidth)**2)
c spinave only (color factors cancel)
      fac=spinave*gw**8*0.5d0*mbsq*(s56-4d0*mb**2)/prop
      qqbWH=fac*s(1,4)*s(2,3)
      qbqWH=fac*s(2,4)*s(1,3)

      do j=-nf,nf
      do k=-nf,nf
      if ((j .gt. 0) .and. (k .lt. 0)) msq(j,k)=Vsq(j,k)*qqbWH
      if ((j .lt. 0) .and. (k .gt. 0)) msq(j,k)=Vsq(j,k)*qbqWH
      enddo
      enddo
      return
      end

