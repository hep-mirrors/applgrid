      subroutine qqb_zz(p,msq)
      implicit none

C----Author R.K.Ellis December 1998
C----modified by JMC to include supplementary diagrams February 1999
c----Matrix element for ZZ production
c----NB: we also include virtual photons
c    in the notation of DKS
C    averaged over initial colours and spins
c    u(-p1)+dbar(-p2)-->e^-(p3) + e^+(p4)   + \mu^-(p5)+ \mu^+(p6)
c    q(-p1)+qbar(-p2)-->l'(p3)  + lbar'(p4) + l(p5)    + lbar(p6)
c    with Z-leptons couplings l1 for (5,6) and l2 for(3,4)
c          and lepton charges q2 for (5,6) and q1 for (3,4)
c    Notation to allow room for p3 --- gluon emission.
c----No statistical factor of 1/2 included.

      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),qdks(mxpart,4),
     . ave,v1(2),v2(2)
     
      double complex qqb(2,2,2),qbq(2,2,2),q_qb,qb_q
      double complex qqb1(2,2,2),qbq1(2,2,2),qqb2(2,2,2),qbq2(2,2,2)
      double complex propz1,propz2,props,a6trees,cprop
      double complex prop12,prop34,prop56

      logical dronly
      double precision FAC
      integer j,k,polq,pol1,pol2
      parameter(ave=0.25d0/xn)

c--- THIS SWITCHES OFF SINGLY RESONANT TERMS REGARDLESS OF ZEROWIDTH
c--- SINCE THEY'RE NOT INCLUDED IN THE REAL EMISSION CONTRIBUTION
      dronly=.true.

      fac=-4D0*esq**2

      v1(1)=l1
      v1(2)=r1
      v2(1)=l2
      v2(2)=r2

c--set msq=0 to initalize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C----Change the momenta to DKS notation 
c   We have --- q(-p1)+qbar(-p2)-->l(p3)+lbar(p4) + l'(p5)+lbar'(p6)
c   DKS have--- q(q2) +qbar(q1) -->mu^-(q3)+mu^+(q4)+e^-(q6)+e^+(q5)

      do j=1,4
      qdks(1,j)=p(1,j)
      qdks(2,j)=p(2,j)
      qdks(3,j)=p(3,j)
      qdks(4,j)=p(4,j)
      qdks(5,j)=p(6,j)
      qdks(6,j)=p(5,j)
      enddo

      call spinoru(6,qdks,za,zb)

c--   s returned from sprod (common block) is 2*dot product

c--   calculate propagators
      if     (dronly  .eqv. .true.) then
      prop12=s(1,2)/dcmplx(s(1,2)-zmass**2,zmass*zwidth)
      prop34=s(3,4)/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=s(5,6)/dcmplx(s(5,6)-zmass**2,zmass*zwidth)
      cprop=dcmplx(1d0)
      elseif (dronly .neqv. .true.) then
      prop12=dcmplx(s(1,2)/(s(1,2)-zmass**2))
      prop34=dcmplx(s(3,4)/(s(3,4)-zmass**2))
      prop56=dcmplx(s(5,6)/(s(5,6)-zmass**2))
      propz1=(s(3,4)-zmass**2)/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      propz2=(s(5,6)-zmass**2)/dcmplx(s(5,6)-zmass**2,zmass*zwidth)
      props=(s(1,2)-zmass**2)/dcmplx(s(1,2)-zmass**2,zmass*zwidth)
      cprop=propz1*propz2*props
      endif
            
c-- here the labels correspond to the polarizations of the
c-- quark, lepton 4 and lepton 6 respectively

      qbq(1,1,1)=A6trees(1,2,6,5,4,3,za,zb) 
      qbq(1,1,2)=A6trees(1,2,6,5,3,4,za,zb) 
      qbq(1,2,1)=A6trees(1,2,5,6,4,3,za,zb) 
      qbq(1,2,2)=A6trees(1,2,5,6,3,4,za,zb) 

      qqb(1,1,1)=A6trees(2,1,6,5,4,3,za,zb)
      qqb(1,1,2)=A6trees(2,1,6,5,3,4,za,zb) 
      qqb(1,2,1)=A6trees(2,1,5,6,4,3,za,zb) 
      qqb(1,2,2)=A6trees(2,1,5,6,3,4,za,zb)
      
      if (dronly .neqv. .true.) then
c---for supplementary diagrams.
      qbq1(1,1,1)=+A6trees(3,4,1,2,5,6,za,zb)
      qbq2(1,1,1)=+A6trees(6,5,1,2,4,3,za,zb)
      qbq1(1,1,2)=-A6trees(4,3,1,2,5,6,za,zb)
      qbq2(1,1,2)=+A6trees(6,5,1,2,3,4,za,zb)      
      qbq1(1,2,1)=+A6trees(3,4,1,2,6,5,za,zb)
      qbq2(1,2,1)=-A6trees(5,6,1,2,4,3,za,zb)
      qbq1(1,2,2)=-A6trees(4,3,1,2,6,5,za,zb)
      qbq2(1,2,2)=-A6trees(5,6,1,2,3,4,za,zb)

      qqb1(1,1,1)=-A6trees(3,4,2,1,5,6,za,zb)
      qqb2(1,1,1)=-A6trees(6,5,2,1,4,3,za,zb)
      qqb1(1,1,2)=+A6trees(4,3,2,1,5,6,za,zb)
      qqb2(1,1,2)=-A6trees(6,5,2,1,3,4,za,zb)      
      qqb1(1,2,1)=-A6trees(3,4,2,1,6,5,za,zb)
      qqb2(1,2,1)=+A6trees(5,6,2,1,4,3,za,zb)
      qqb1(1,2,2)=+A6trees(4,3,2,1,6,5,za,zb)
      qqb2(1,2,2)=+A6trees(5,6,2,1,3,4,za,zb)
      endif

      do j=1,2
      do k=1,2
      qbq(2,j,k)=-qqb(1,j,k)
      qqb(2,j,k)=-qbq(1,j,k)
      qbq1(2,j,k)=-qqb1(1,j,k)
      qqb1(2,j,k)=-qbq1(1,j,k)
      qbq2(2,j,k)=-qqb2(1,j,k)
      qqb2(2,j,k)=-qbq2(1,j,k)
      enddo
      enddo

      do j=-nf,nf
      k=-j
      msq(j,k)=0d0

      if (j.eq.0) go to 20

c      if (abs(j) .le. 2) write(*,*) 'MCFM, j=',j

      if ((j .gt. 0).and.(k .lt. 0)) then
      do polq=1,2
      do pol1=1,2
      do pol2=1,2
      if     (polq .eq. 1) then
       q_qb=(prop56*v2(pol1)*l(j)+q2*q(j))
     .     *(prop34*v1(pol2)*l(j)+q1*q(j))*qqb(polq,pol1,pol2)
        if (dronly .neqv. .true.) then
         q_qb=q_qb
     .       +(prop56*v1(pol2)*v2(pol1)+q1*q2)
     .       *(prop12+v1(pol2)*l(j)+q1*q(j))*qqb1(polq,pol1,pol2)
     .       +(prop34*v1(pol2)*v2(pol1)+q1*q2)
     .       *(prop12*v2(pol1)*l(j)+q2*q(j))*qqb2(polq,pol1,pol2)
        endif
      elseif (polq .eq. 2) then
       q_qb=(prop56*v2(pol1)*r(j)+q2*q(j))
     .     *(prop34*v1(pol2)*r(j)+q1*q(j))*qqb(polq,pol1,pol2)
        if (dronly .neqv. .true.) then
         q_qb=q_qb
     .       +(prop56*v1(pol2)*v2(pol1)+q1*q2)
     .       *(prop12+v1(pol2)*r(j)+q1*q(j))*qqb1(polq,pol1,pol2)
     .       +(prop34*v1(pol2)*v2(pol1)+q1*q2)
     .       *(prop12*v2(pol1)*r(j)+q2*q(j))*qqb2(polq,pol1,pol2)
        endif
      endif
C-- Inclusion of width a la Baur and Zeppenfeld
      q_qb=FAC*q_qb*cprop
      msq(j,k)=msq(j,k)+ave*abs(q_qb)**2

c      if (abs(j) .le. 2) 
c     .write(*,*)'(',polq,', ',pol1,', ',pol2,') = ',4d0*ave*abs(q_qb)**2

      enddo
      enddo
      enddo

      elseif ((j .lt. 0).and.(k .gt. 0)) then

      do polq=1,2
      do pol1=1,2
      do pol2=1,2
      if     (polq .eq. 1) then
       qb_q=(prop56*v2(pol1)*l(k)+q2*q(k))
     .     *(prop34*v1(pol2)*l(k)+q1*q(k))*qbq(polq,pol1,pol2)
        if (dronly .neqv. .true.) then
         qb_q=qb_q
     .       +(prop56*v1(pol2)*v2(pol1)+q1*q2)
     .       *(prop12+v1(pol2)*l(k)+q1*q(k))*qbq1(polq,pol1,pol2)
     .       +(prop34*v1(pol2)*v2(pol1)+q1*q2)
     .       *(prop12*v2(pol1)*l(k)+q2*q(k))*qbq2(polq,pol1,pol2)
        endif
      elseif (polq .eq. 2) then
       qb_q=(prop56*v2(pol1)*r(k)+q2*q(k))
     .     *(prop34*v1(pol2)*r(k)+q1*q(k))*qbq(polq,pol1,pol2)
        if (dronly .neqv. .true.) then
         qb_q=qb_q
     .       +(prop56*v1(pol2)*v2(pol1)+q1*q2)
     .       *(prop12+v1(pol2)*r(k)+q1*q(k))*qbq1(polq,pol1,pol2)
     .       +(prop34*v1(pol2)*v2(pol1)+q1*q2)
     .       *(prop12*v2(pol1)*r(k)+q2*q(k))*qbq2(polq,pol1,pol2)
        endif
      endif
C-- Inclusion of width a la Baur and Zeppenfeld
      qb_q=FAC*qb_q*cprop
      msq(j,k)=msq(j,k)+ave*abs(qb_q)**2

c      if (abs(j) .le. 2) 
c     .write(*,*)'(',polq,', ',pol1,', ',pol2,') = ',4d0*ave*abs(qb_q)**2

      enddo
      enddo
      enddo

      endif

 20   continue
      enddo

      return
      end
