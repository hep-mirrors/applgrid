      subroutine dips_mass(nd,p,ip,jp,kp,sub,subv,msq,msqv,
     . subr_born,subr_corr)
      implicit none
************************************************************************
*     Author: Keith Ellis                                              *
*     June 2002                                                        *
*     Calculates the nj-jet subtraction term corresponding to dipole   *
*     nd with momentum p and dipole kinematics (ip,jp) wrt kp          *
*     Automatically chooses dipole kind                                *
*     Returns the dipoles in sub,subv and matrix elements in msq,msqv  *
*     nd labels the dipole configurations                              *
*     ip labels the emitter parton                                     *
*     jp labels the emitted parton                                     *
*     kp labels the spectator parton                                   *
*     subr_born is the subroutine which call the born process          *
*     subr_corr is the subroutine which call the born process dotted   *
*     with vec for an emitted gluon only                               *
************************************************************************
      include 'constants.f'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'process.f'
      include 'alfacut.f'
      double precision p(mxpart,4),ptrans(mxpart,4),sub(4),subv,vecsq,
     . x,omx,z,omz,y,omy,u,omu,pij,pik,pjk,dot,q(4),qsq,qij(4),qijsq,
     . vec(4),root,vtilde,pold(mxpart,4),pext(mxpart,4)
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),zp,zm
      double precision mksq,misq,mjsq,mijsq,muisq,mujsq,muksq,
     . muijsq,kappa,vijk,vtijk,viji,ztmi,ztmj,muk,mqsq
      double precision yp,mass2,width2,mass3,width3
      integer nd,ip,jp,kp,nu,j,jproc,n2,n3
      logical incldip(0:maxd)
      common/incldip/incldip
      common/breit/n2,n3,mass2,width2,mass3,width3
      external subr_born,subr_corr
      parameter(kappa=0d0)
      logical first
      data first/.true./

      if (first) then
        first=.false.
        write(6,*) 'dips_mass:mass2',mass2
      endif 
      mqsq=mass2**2

      do nu=1,4
      do j=1,mxpart
      ptrans(j,nu)=0d0
      pold(j,nu)=p(j,nu)
      enddo
      enddo

      if ((case .eq. 't_bbar') .or. (case .eq. 'bq_tpq')) then
c--- if we're doing single-top, reduce # of momenta from 7 to 5 
        do nu=1,4
          p(3,nu)=pold(3,nu)+pold(4,nu)+pold(5,nu)
          p(4,nu)=pold(6,nu)
          p(5,nu)=pold(7,nu)
          p(6,nu)=0d0
          p(7,nu)=0d0
        enddo
      endif
      
      if ((case .eq. 'W_twdk') .or. (case .eq. 'W_cwdk')) then
c--- if we're doing W+t, reduce # of momenta from 8 to 6 
        do nu=1,4
          p(3,nu)=pold(3,nu)
          p(4,nu)=pold(4,nu)
          p(5,nu)=pold(5,nu)+pold(6,nu)+pold(7,nu)
          p(6,nu)=pold(8,nu)
          p(7,nu)=0d0
          p(8,nu)=0d0
        enddo
      endif
      
C---Initialize the dipoles to zero
      do j=1,4
      sub(j)=0d0
      enddo
      subv=0d0
      call zeromsq(msq,msqv)
      incldip(nd)=.true.

C--- default is all particles massless
      misq=0d0
      mjsq=0d0
      mksq=0d0
      mijsq=0d0

      pij=two*dot(p,ip,jp)
      pik=two*dot(p,ip,kp)
      pjk=two*dot(p,jp,kp)

      if ((ip .le. 2) .and. (kp .le. 2)) then
***********************************************************************
*************************** INITIAL-INITIAL ***************************
***********************************************************************
        omx=-(pij+pjk)/pik
        x=one-omx
        vtilde=pij/pik
                
C---Modification so that only close to singular subtracted
        if (-vtilde .gt. aii) then
           incldip(nd)=.false.
           goto 99
        endif
        
        call transform_mass(p,ptrans,x,ip,jp,kp,misq,mjsq,mksq,mijsq)

        if ((case .eq. 't_bbar') .or. (case .eq. 'bq_tpq')
     .  .or.(case .eq. 'W_twdk') .or. (case .eq. 'W_cwdk')) then
          if ((case .eq. 'W_twdk') .or. (case .eq. 'W_cwdk')) then
            call extend_trans_wt(pold,p,ptrans,pext)
          else
            call extend_trans(pold,p,ptrans,pext)
          endif
          do j=1,mxpart
          do nu=1,4
            ptrans(j,nu)=pext(j,nu)
          enddo
          enddo
        endif

        call storeptilde(nd,ptrans)
        
        do nu=1,4
          vec(nu)=p(jp,nu)-pij/pik*p(kp,nu)
        enddo
        vecsq=-pij*pjk/pik
        call subr_born(ptrans,msq)
        call subr_corr(ptrans,vec,ip,msqv)

        sub(qq)=-gsq/x/pij*(two/omx-one-x)
        sub(gq)=-gsq/pij
        sub(qg)=-gsq/x/pij*(one-two*x*omx)
        sub(gg)=-2d0*gsq/x/pij*(x/omx+x*omx)
        subv   =+4d0*gsq/x/pij*omx/x/vecsq

***********************************************************************
*************************** INITIAL-FINAL *****************************
***********************************************************************
      elseif ((ip .le. 2) .and. (kp .gt. 2)) then
        omx=-pjk/(pij+pik)
        x=one-omx
        u=pij/(pij+pik)
        omu=pik/(pij+pik)

C---determine mass of spectator
c      mksq=max(p(kp,4)**2-p(kp,1)**2-p(kp,2)**2-p(kp,3)**2,0d0)
c      if (mksq.gt.0d0) then
c         muksq=mksq/2d0/(-dot(p,ip,jp)-dot(p,ip,kp))
c         zp=omx/(omx+muksq)
c      else
c         zp=1d0
c      endif     
      
C---Modification so that only close to singular subtracted
        if (u .gt. aif) goto 99
        
C---npart is the number of particles in the final state
C---transform the momenta so that only the first npart+1 are filled
        call transform_mass(p,ptrans,x,ip,jp,kp,misq,mjsq,mksq,mijsq)

        if ((case .eq. 't_bbar') .or. (case .eq. 'bq_tpq')
     .  .or.(case .eq. 'W_twdk') .or. (case .eq. 'W_cwdk')) then
          if ((case .eq. 'W_twdk') .or. (case .eq. 'W_cwdk')) then
            call extend_trans_wt(pold,p,ptrans,pext)
          else
            call extend_trans(pold,p,ptrans,pext)
          endif
          do j=1,mxpart
          do nu=1,4
            ptrans(j,nu)=pext(j,nu)
          enddo
          enddo
        endif

        call storeptilde(nd,ptrans)
        do nu=1,4
           vec(nu)=p(jp,nu)/u-p(kp,nu)/omu
        enddo
        call subr_born(ptrans,msq)
        call subr_corr(ptrans,vec,ip,msqv)        
        sub(qq)=-gsq/x/pij*(two/(omx+u)-one-x)
        sub(gq)=-gsq/pij
        sub(qg)=-gsq/x/pij*(one-two*x*omx)
        sub(gg)=-2d0*gsq/x/pij*(one/(omx+u)-one+x*omx)
        subv   =-4d0*gsq/x/pij*(omx/x*u*(one-u)/pjk)
***********************************************************************
*************************** FINAL-INITIAL *****************************
***********************************************************************
      elseif ((ip .gt. 2) .and. (kp .le. 2)) then
        do jproc=1,4
        if ((jproc.eq.qq) .and. (qqproc .eqv. .false.)) goto 79
        if ((jproc.eq.gq) .and. (gqproc .eqv. .false.)) goto 79
        if ((jproc.eq.qg) .and. (qgproc .eqv. .false.)) goto 79
        if ((jproc.eq.gg) .and. (ggproc .eqv. .false.)) goto 79

        if (jproc.eq.qq) then
        mijsq=mqsq
c--- the masses of i and j have been switched
        misq=mqsq
        mjsq=0d0
        elseif (jproc.eq.qg) then
        go to 79
        elseif (jproc.eq.gq) then
        mijsq=0d0
        misq=mqsq
        mjsq=mqsq
        elseif (jproc.eq.gg) then
        goto 79
        endif
        omx=(mijsq-misq-mjsq-pij)/(pjk+pik)
        x=one-omx

c---Modification so that only close to singular subtracted
        if (omx .gt. afi) goto 99

        do nu=1,4
          qij(nu)=p(ip,nu)+p(jp,nu)
          q(nu)=qij(nu)+p(kp,nu)
        enddo
        qsq=q(4)**2-q(1)**2-q(2)**2-q(3)**2
        qijsq=qij(4)**2-qij(1)**2-qij(2)**2-qij(3)**2
        call transform_mass(p,ptrans,x,ip,jp,kp,misq,mjsq,mksq,mijsq)

        if ((case .eq. 't_bbar') .or. (case .eq. 'bq_tpq')
     .  .or.(case .eq. 'W_twdk') .or. (case .eq. 'W_cwdk')) then
          if ((case .eq. 'W_twdk') .or. (case .eq. 'W_cwdk')) then
            call extend_trans_wt(pold,p,ptrans,pext)
          else
            call extend_trans(pold,p,ptrans,pext)
          endif
          do j=1,mxpart
          do nu=1,4
            ptrans(j,nu)=pext(j,nu)
          enddo
          enddo
        endif

        call storeptilde(nd,ptrans)
        z=pik/(pik+pjk)
        omz=pjk/(pik+pjk)
c--- note that musq is related to msq by musq = msq/(2pij_tilde.pa)
c--- and 2pij_tilde.pa = (Qsq-mijsq)/x
        zp=omx*(Qsq-mijsq)/x+mijsq+misq-mjsq
        root=dsqrt(zp**2-4d0*misq*mjsq)
        zm=(zp-root)/(2d0*(omx*(Qsq-mijsq)/x+mijsq))
        zp=(zp+root)/(2d0*(omx*(Qsq-mijsq)/x+mijsq))

        call subr_born(ptrans,msq)

c        do nu=1,4
c          vec(nu)=z*p(ip,nu)-omz*p(jp,nu)
c        enddo

c        call subr_corr(ptrans,vec,5,msqv)

        if (jproc .eq. qq) then
        sub(qq)=+gsq/x/(qijsq-mijsq)*(two/(omz+omx)-one-z-2d0*mqsq/pij)
        elseif (jproc .eq. gq) then
        sub(gq)=+gsq/x/(qijsq-mijsq)
        subv   =+4d0*gsq/x/qijsq/(qijsq-mijsq)
        endif
 79     continue
        enddo

***********************************************************************
**************************** FINAL-FINAL ******************************
***********************************************************************
      elseif ((ip .gt. 2) .and. (kp .gt. 2)) then
c------Eq-(5.2)    

C----Form momentum vectors
      y=pij/(pij+pjk+pik)
      z=pik/(pjk+pik)
      omz=one-z
      omy=one-y
       do nu=1,4
C  create 4-momentum of sum emitter+emittee
          qij(nu)=p(ip,nu)+p(jp,nu)
C  create 4-momentum of sum emitter+emittee+spectator 
          q(nu)=qij(nu)+p(kp,nu)
       enddo
C  and square them 
      Qsq=q(4)**2-q(1)**2-q(2)**2-q(3)**2
      qijsq=qij(4)**2-qij(1)**2-qij(2)**2-qij(3)**2



C---determine mass of spectator
      mksq=max(p(kp,4)**2-p(kp,1)**2-p(kp,2)**2-p(kp,3)**2,0d0)
      if (mksq.gt.0d0) then
         muk=sqrt(mksq/qsq)
         yp=(1d0-muk)/(1d0+muk)
      else
         yp=1d0
      endif     

      if (y .gt. yp) then
         write(6,*) 'Problems with phase space in dips_mass.f'
         stop
      endif

C---Modification so that only close to singular subtracted
       if (y .gt. aff*yp) then
         incldip(nd)=.false.
         go to 99
       endif

C---loop over the different possibilities which have different kinematics
      do jproc=1,4
      if ((jproc.eq.qq) .and. (qqproc .eqv. .false.)) goto 80
      if ((jproc.eq.gq) .and. (gqproc .eqv. .false.)) goto 80
      if ((jproc.eq.qg) .and. (qgproc .eqv. .false.)) goto 80
      if ((jproc.eq.gg) .and. (ggproc .eqv. .false.)) goto 80


      if (jproc.eq.qq) then
C q->qg
      mijsq=mqsq
c--- the masses of i and j have been switched
      misq=mqsq
      mjsq=0d0
      elseif (jproc.eq.qg) then
C q->gq
      go to 80
      elseif (jproc.eq.gq) then
C g->qqbar
      mijsq=0d0
      misq=mqsq
      mjsq=mqsq
      elseif (jproc.eq.gg) then
C g->gg
      mijsq=0d0
      misq=0d0
      mjsq=0d0
      endif
      

      muisq=misq/Qsq
      mujsq=mjsq/Qsq
      muksq=mksq/Qsq
      muijsq=mijsq/Qsq
      muk=dsqrt(muksq)


c      viji=sqrt((1d0-muijsq-muisq)**2-4d0*mijsq*muisq)
c     . /(1d0-muijsq-muisq)
c      write(6,*) 'viji',viji
c      vijk=sqrt((one-qijsq/Qsq-muksq)**2-4d0*qijsq/Qsq*muksq)
c     . /(one-qijsq/Qsq-muksq)
c      write(6,*) vijk

      viji=dsqrt(((1d0-mujsq-muisq-muksq)*y)**2-4d0*muisq*mujsq)
     . /((1d0-mujsq-muisq-muksq)*y+2d0*muisq)

      vijk=dsqrt((2d0*muksq+(1d0-mujsq-muisq-muksq)*omy)**2-4d0*muksq)
     . /((1d0-mujsq-muisq-muksq)*omy)


c      ym=2d0*mui*muj/(1d0-muisq-mujsq-muksq)
      yp=1d0-2d0*muk*(1d0-muk)/(1d0-muisq-mujsq-muksq)

      zp=(2d0*muisq+(1d0-muisq-mujsq-muksq)*y)
     . /(2d0*(muisq+mujsq+(1d0-muisq-mujsq-muksq)*y)) 
      zm=zp*(1d0-viji*vijk)
      zp=zp*(1d0+viji*vijk)
C---calculate the ptrans-momenta 
       call transform_mass(p,ptrans,y,ip,jp,kp,misq,mjsq,mksq,mijsq)

        if ((case .eq. 't_bbar') .or. (case .eq. 'bq_tpq')
     .  .or.(case .eq. 'W_twdk') .or. (case .eq. 'W_cwdk')) then
          if ((case .eq. 'W_twdk') .or. (case .eq. 'W_cwdk')) then
            call extend_trans_wt(pold,p,ptrans,pext)
          else
            call extend_trans(pold,p,ptrans,pext)
          endif
          do j=1,mxpart
          do nu=1,4
            ptrans(j,nu)=pext(j,nu)
          enddo
          enddo
        endif

c       write(6,*) 'Dipole ',nd, 'ptrans'
c       call writeout(ptrans)

C have to enhance the store so that it works
       call storeptilde(nd,ptrans)

       ztmi=z-0.5d0+0.5d0*vijk
       ztmj=omz-0.5d0+0.5d0*vijk

       call subr_born(ptrans,msq)

       do nu=1,4
         vec(nu)=ztmi*p(ip,nu)-ztmj*p(jp,nu)
       enddo

       if (ip .lt. kp) then
         call subr_corr(ptrans,vec,5,msqv)
       else
         call subr_corr(ptrans,vec,6,msqv)
       endif
              
      if (jproc .eq. qq) then
        vtijk=sqrt((one-muijsq-muksq)**2-4d0*muijsq*muksq)
     .  /(one-muijsq-muksq)

        sub(qq)=gsq/(qijsq-mijsq)*(two/(one-z*omy)
     .  -vtijk/vijk*(one+z+2d0*mqsq/pij))
     
c        write(6,*) 'velocities ratio, vtijk/vijk=',vtijk/vijk
c        write(6,*) 'sub(qq)=',sub(qq)
C---debug
c       sub(qq)=gsq/pij*(two/(one-omz*omy)
c     . -(one+omz+2d0*mqsq/pij))
C---debug
      elseif (jproc .eq. gq) then
        sub(gq)=gsq/(qijsq-mijsq)*(one-two*kappa*(zp*zm-mqsq/qijsq))
        subv   =+4d0*gsq/pij/pij
      elseif (jproc .eq. gg) then
        vijk=sqrt((one-qijsq/Qsq-muksq)**2-4d0*qijsq/Qsq*muksq)
     .  /(one-qijsq/Qsq-muksq)
        sub(gg)=two*gsq/(qijsq-mijsq)*(one/(one-z*omy)+one/(one-omz*omy)
     .    -(two-kappa*zp*zm)/vijk)
      subv   =+4d0*gsq/(qijsq-mijsq)/pij/vijk
      endif
 80   continue
      enddo
      endif
      
c--- fall through to here, so that p retains the value it entered with
   99 continue   
      
      do j=1,mxpart
      do nu=1,4
        p(j,nu)=pold(j,nu)
      enddo
      enddo

      return
      end
      
