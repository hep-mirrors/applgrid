      subroutine gg_hg_gvec(p,n,in,msq)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
C  in is the label of the momentum contracted with n
      integer j,k,in
      double precision msq(-nf:nf,-nf:nf)
      double precision n(4),p(mxpart,4),dot,hdecay,s34,fac,
     . qqghn,ggghn,p1p2(-1:1,-1:1)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C   Deal with Higgs decay to b-bbar
      s34=2d0*Dot(p,3,4)+2d0*mb**2
      hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s34-4d0*mb**2) 
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      fac=hdecay

      do j=-1,+1
      do k=-1,+1
      p1p2(j,k)=0d0
      enddo
      enddo

      if (in .eq. 1) then
      p1p2(0,-1)=-aveqg*fac*qqghn(2,5,1,p,n)
      p1p2(0,+1)=-aveqg*fac*qqghn(2,5,1,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(5,2,1,p,n)
      elseif (in .eq. 2) then
      p1p2(+1,0)=-aveqg*fac*qqghn(1,5,2,p,n)
      p1p2(-1,0)=-aveqg*fac*qqghn(5,1,2,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(1,5,2,p,n)
      elseif (in .eq. 5) then     
      p1p2(1,-1)=+aveqq*fac*qqghn(1,2,5,p,n)
      p1p2(-1,1)=+aveqq*fac*qqghn(2,1,5,p,n)
      p1p2(0,0)=+avegg*fac*ggghn(1,2,5,p,n)
      endif

      do j=-nf,nf
      do k=-nf,nf
      if     ((j .gt. 0) .and. (k .eq. -j)) then
          msq(j,k)=p1p2(1,-1)
      elseif ((j .lt. 0) .and. (k .eq. -j)) then
          msq(j,k)=p1p2(-1,1)
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
          msq(j,k)=p1p2(0,0)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    p1p2(+1,0)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    p1p2(-1,0)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=
     &    p1p2(0,+1)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=
     &    p1p2(0,-1)
      endif
      enddo
      enddo
 
      return
      end

      double precision function qqghn(j1,j2,j5,p,n)
      implicit none
C---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> H((p3+p4)+g(p5)
c   contracted with the vector n(mu)
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j1,j2,j5
      double precision Asq,p(mxpart,4),n(4),nDn,nDp1,nDp2,nDp5,
     . dot,s,t,u

      nDp1=n(4)*p(j1,4)-n(3)*p(j1,3)-n(2)*p(j1,2)-n(1)*p(j1,1)
      nDp2=n(4)*p(j2,4)-n(3)*p(j2,3)-n(2)*p(j2,2)-n(1)*p(j2,1)
      nDp5=n(4)*p(j5,4)-n(3)*p(j5,3)-n(2)*p(j5,2)-n(1)*p(j5,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

c--- appropriate scale is approx 1d-3*energy(incoming)
c--- so of order(1) for the Tevatron
      if (abs(nDp5).gt.1d-3*abs(p(1,4))) then 
         write(*,*) 'Error for :',j1,j2,j5
         write(*,*) 'cutoff',1d-3*abs(p(1,4))
         write(6,*) 'nDp5',nDp5
         call flush(6)
         stop
      endif

      Asq=(as/(3d0*pi))**2/vevsq

      s=2d0*Dot(p,j1,j2)
      t=2d0*Dot(p,j1,j5)
      u=2d0*Dot(p,j2,j5)

c--- RKE Schoonship answer, ncalc.out
c     Id,bit=g^2*A^2*V/2*(2*(nDp1*u-nDp2*t)^2/s^2
c           + 0.5*nDn*(u+t)^2/s)

      qqghn=-Asq*gsq*V/2d0*(2d0*(nDp1*u-nDp2*t)**2/s**2
     .                    +0.5d0*nDn*(u+t)**2/s)

      return
      end

      double precision function ggghn(j1,j2,j5,p,n)
      implicit none
C---calculates the amplitude squared for the process
c   g(p1)+g(p2) --> H((p3+p4)+g(p5)
c   contracted with the vector n(mu)
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j1,j2,j5
      double precision Asq,p(mxpart,4),n(4),nDn,nDp1,nDp2,nDp5,
     . dot,s,t,u,sh

      nDp1=n(4)*p(j1,4)-n(3)*p(j1,3)-n(2)*p(j1,2)-n(1)*p(j1,1)
      nDp2=n(4)*p(j2,4)-n(3)*p(j2,3)-n(2)*p(j2,2)-n(1)*p(j2,1)
      nDp5=n(4)*p(j5,4)-n(3)*p(j5,3)-n(2)*p(j5,2)-n(1)*p(j5,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

c--- appropriate scale is approx 1d-3*energy(incoming)
c--- so of order(1) for the Tevatron
      if (abs(nDp5).gt.1d-3*abs(p(1,4))) then 
         write(*,*) 'Error for :',j1,j2,j5
         write(*,*) 'cutoff',1d-3*abs(p(1,4))
         write(6,*) 'nDp5',nDp5
         call flush(6)
         stop
      endif

      Asq=(as/(3d0*pi))**2/vevsq

      s=2d0*Dot(p,j1,j2)
      t=2d0*Dot(p,j1,j5)
      u=2d0*Dot(p,j2,j5)
      sh=s+t+u
      
c--- JMC answer, gggH.frm
c -f(b,a,c)^2/s12/s13/s23*(
c -(n.n)/2*(s12^4+s13^4+s23^4+mHsq^4-2*(s13^2*s23^2+s12^2*mHsq^2))
c +2*(p1.n*s23-p2.n*s13)^2/s13/s23*(s13^2*s23^2+s12^2*mHsq^2)/s12);

      ggghn=Asq*gsq*V*xn*(
     . -nDn/2d0*(s**4+t**4+u**4+sh**4-2d0*(t**2*u**2+s**2*sh**2))
     . +2d0*(nDp1*u-nDp2*t)**2/t/u*(t**2*u**2+s**2*sh**2)/s)/s/t/u

      return
      end
      
