      subroutine w2jetnx(i1,i2,i3,i4,i5,i6,p,n,za,zb,zab,zba)
C----matrix element squared with p5 line contracted with n(mu)
C----nDp6 should be equal to zero
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'mmsqv_cs.f'
      double complex qcdabn(2,2,2),qcdban(2,2,2),qedn(2,2,2)
      double complex zab(mxpart,mxpart),zba(mxpart,mxpart)
      double precision msq1,msq2,msqq,n(4),p(mxpart,4)
      double precision nDp5,nDp6
      integer i1,i2,i3,i4,i5,i6

      nDp5=n(4)*p(i5,4)-n(3)*p(i5,3)-n(2)*p(i5,2)-n(1)*p(i5,1)
      nDp6=n(4)*p(i6,4)-n(3)*p(i6,3)-n(2)*p(i6,2)-n(1)*p(i6,1)
c--- appropriate scale is approx 1d-3*energy(incoming)
c--- so of order(1) for the Tevatron
      if (abs(nDp6).gt.1d-3*abs(p(i1,4))) then 
         write(*,*) 'Error for :',i1,i2,i3,i4,i5,i6
         write(*,*) 'cutoff',1d-3*abs(p(i1,4))
         write(6,*) 'nDp5',nDp5
         write(6,*) 'nDp6',nDp6
         call flush(6)
         stop
      endif

      call subqcdn(i1,i2,i3,i4,i5,i6,nDp5,za,zb,zab,zba,qcdabn,qcdban)
            
C--first argument is quark line
C--second argument is polarization of i5 line
C  1=L,2=R
      qedn(1,1,1)=qcdabn(1,1,1)+qcdban(1,1,1) 
      qedn(2,1,1)=qcdabn(2,1,1)+qcdban(2,1,1) 

      msq1= abs(qcdabn(1,1,1))**2+abs(qcdabn(2,1,1))**2
      msq2= abs(qcdban(1,1,1))**2+abs(qcdban(2,1,1))**2 
      msqq= abs(qedn(1,1,1))**2+abs(qedn(2,1,1))**2

      mmsqv_cs(0,+1,+1)=-ninth*msqq
      mmsqv_cs(1,+1,+1)=msq1
      mmsqv_cs(2,+1,+1)=msq2

      return
      end

