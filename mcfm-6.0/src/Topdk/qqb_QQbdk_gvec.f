      subroutine qqb_QQbdk_gvec(p,n,in,msqv)
      implicit none

c----Matrix element for tt production
C----averaged over initial colours and spins
c    line in contracted with the vector n(mu)
c     g(-p1)+g(-p2)--> t(p3,p4,p5)+tb(p6,p7,p8)

      include 'constants.f'
      include 'qcdcouple.f'
      include 'msqv_cs.f'

C in is the label of the contracted line
      integer j,k,in,icol,nu
      double precision q1(4),q2(4),q3(4),q4(4),q5(4),q6(4),q7(4),q8(4)
      double precision msqv(-nf:nf,-nf:nf),p(mxpart,4),res(0:2)
      double precision n(4),nDp1

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0d0
      do icol=0,2
        msqv_cs(icol,j,k)=0d0
      enddo
      enddo
      enddo

      nDp1=n(4)*p(in,4)-n(1)*p(in,1)-n(2)*p(in,2)-n(3)*p(in,3)
      if (abs(nDp1).gt.1d-3*abs(p(1,4))) then 
         write(*,*) 'Error for in=:',in
         write(*,*) 'cutoff',1d-3*abs(p(1,4))
         write(6,*) 'qqb_QQbdk_gvec:nDp1',nDp1
         call flush(6)
         stop
      endif


      do nu=1,4
      q1(nu)=p(1,nu)
      q2(nu)=p(2,nu)
      q3(nu)=p(3,nu)
      q4(nu)=p(4,nu)
      q5(nu)=p(5,nu)
      q6(nu)=p(6,nu)
      q7(nu)=p(7,nu)
      q8(nu)=p(8,nu)
      enddo

c      call ampsqggQQbdkn(n,in,q1,q2,q3,q4,q5,q6,q7,q8,res)
      call KMampsqggQQbdkn(n,in,q1,q2,q3,q4,q5,q6,q7,q8,res)

      do icol=0,2
      msqv_cs(icol,0,0)=gsq**2*avegg*res(icol)
      enddo

      msqv(0,0)=msqv_cs(0,0,0)+msqv_cs(1,0,0)+msqv_cs(2,0,0)


      return
      end
      
