      subroutine writeinfo(unitno,xsec,xsec_err,itno)
************************************************************************
*   Routine to write out run information to a desired unit             *
************************************************************************
      implicit none
      include 'constants.f'
      include 'removebr.f'
      include 'PDFerrors.f'
      include 'process.f'
      integer unitno,j,k,itno
      double precision xsec,xsec_err
      double precision lord_bypart(-1:1,-1:1),lordnorm,rescale
      double precision ggpart,gqpart,qgpart,qqpart,qqbpart,
     . gqbpart,qbgpart,qbqbpart,qbqpart
      
      character*4 part
      character*30 runstring
      logical creatent,dswhisto,dryrun,makecuts
      integer nproc,ih1,ih2,itmx1,itmx2,ncall1,ncall2,origij
      integer NPTYPE,NGROUP,NSET
      double precision sqrts
      double precision Rcut
 
      common/outputflags/creatent,dswhisto      

      common/nproc/nproc
      common/part/part
      common/runstring/runstring
      common/energy/sqrts
      common/density/ih1,ih2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/dryrun/dryrun
      
      common/pdflib/NPTYPE,NGROUP,NSET
      
      common/Rcut/Rcut
      common/makecuts/makecuts

      common/origij/origij

      common/bypart/lord_bypart
      common/finalpart/ggpart,gqpart,qgpart,qqpart,qqbpart

      if (itno .gt. 0) then
c--- write warning that result is only intermediate; populate the
c--- variables in finalpart (normally done in mcfm_exit)
      write(unitno,*) '( Intermediate result for iteration',itno,')'
      lordnorm=0d0
      do j=-1,1
      do k=-1,1
        lordnorm=lordnorm+lord_bypart(j,k)
      enddo
      enddo
      ggpart=lord_bypart( 0, 0)/lordnorm
      gqpart=lord_bypart( 0,+1)/lordnorm
      gqbpart=lord_bypart( 0,-1)/lordnorm
      qgpart=lord_bypart(+1, 0)/lordnorm
      qbgpart=lord_bypart(-1, 0)/lordnorm
      qqpart=lord_bypart(+1,+1)/lordnorm
      qbqbpart=lord_bypart(-1,-1)/lordnorm
      qqbpart=lord_bypart(+1,-1)/lordnorm
      qbqpart=lord_bypart(-1,+1)/lordnorm
      endif
      write(unitno,55) '( Cross-section is: ',xsec,' +/-',xsec_err,')'
      write(unitno,*)

c--- for gg->H+X processes, also write out the cross section
c---  normalized by sigma(gg->H, finite mt)/sigma(gg->H, mt-> infinity)
      if ( (case(1:5) .eq. 'ggfus') .or. (case(1:3) .eq. 'HWW')
     . .or.(case(1:3) .eq. 'HZZ')) then
        call finitemtcorr(rescale)
        write(unitno,55)'( Rescaled x-sec is: ',xsec*rescale,' +/-',
     .                                          xsec_err*rescale,')'
        write(unitno,*)
      endif
     
      write(unitno,*) '( Contribution from parton sub-processes:'
      write(unitno,95) '   GG    ',ggpart*xsec,ggpart*100d0
      write(unitno,95) '   GQ    ',gqpart*xsec,gqpart*100d0
      write(unitno,95) '   GQB   ',gqbpart*xsec,gqbpart*100d0
      write(unitno,95) '   QG    ',qgpart*xsec,qgpart*100d0
      write(unitno,95) '   QBG   ',qbgpart*xsec,qbgpart*100d0
      write(unitno,95) '   QQ    ',qqpart*xsec,qqpart*100d0
      write(unitno,95) '   QBQB  ',qbqbpart*xsec,qbqbpart*100d0
      write(unitno,95) '   QQB   ',qqbpart*xsec,qqbpart*100d0
      write(unitno,95) '   QBQ   ',qbqpart*xsec,qbqpart*100d0
      write(unitno,*)

      if (PDFerrors) then
        do j=0,maxPDFsets
          write(unitno,56) j,PDFxsec(j)
        enddo
        write(unitno,*)
      endif

c--- new routine for writing out contents of input file
      call writeinput(unitno,' (',' )','WRITEALL')


c--- old lines for writing out inputs

c      write(unitno,*) '( Run corresponds to this input file)'
c      write(unitno,*)
c      write(unitno,*)
c     . '( [Flags to specify the mode in which MCFM is run] )'
c      write(unitno,98) evtgen,'evtgen'
c      write(unitno,98) creatent,'creatent'
c      write(unitno,98) skipnt,'skipnt'
c      write(unitno,98) dswhisto,'dswhisto'
c
c      write(unitno,*)
c      write(unitno,*)
c     . '( [General options to specify the process and execution] )'
c      write(unitno,97) nproc,'nproc'
c      write(unitno,96) part,'part'
c      write(unitno,96) runstring,'runstring'
c      write(unitno,99) sqrts,'sqrts'
c      write(unitno,97) ih1,'ih1'
c      write(unitno,97) ih2,'ih2'
c      write(unitno,99) hmass,'hmass'
cc--- catch special scale choices for stop+b process
c      if ( (nproc .eq. 231) .or. (nproc .eq. 236)      
c     . .or.(nproc .eq. 241) .or. (nproc .eq. 246)      
c     . .or.(nproc .eq. 242) .or. (nproc .eq. 247) ) then
c         write(unitno,99) renscale_L,'renscale_L'
c         write(unitno,99) facscale_L,'facscale_L'
c         write(unitno,99) renscale_H,'renscale_H'
c         write(unitno,99) facscale_H,'facscale_H'
c      else
c         write(unitno,99) scale,'scale'
c         write(unitno,99) facscale,'facscale'
c      endif
c
c      write(unitno,98) dynamicscale,'dynamicscale'
c      write(unitno,98) zerowidth,'zerowidth'
c      write(unitno,98) removebr,'removebr'
c      write(unitno,97) itmx1,'itmx1'
c      write(unitno,97) ncall1,'ncall1'
c      write(unitno,97) itmx2,'itmx2'
c      write(unitno,97) ncall2,'ncall2'
c      write(unitno,97) origij,'ij'
c      write(unitno,98) dryrun,'dryrun'
c      write(unitno,98) Qflag,'Qflag'
c      write(unitno,98) Gflag,'Gflag'
c      
c      write(unitno,*)
c      write(unitno,*) 
c     . '( [Heavy quark masses] )'
c      write(unitno,99) mt,'top mass'
c      write(unitno,99) mb,'bottom mass'
c      write(unitno,99) mc,'charm mass'
c
c      write(unitno,*)
c      write(unitno,*) 
c     . '( [Pdf selection] )'
c      write(unitno,96) pdlabel,'pdlabel '
c      write(unitno,97) NGROUP,'NGROUP'
c      write(unitno,97) NSET,'NSET'
c      write(unitno,96) PDFname,'LHAPDF group'
c      write(unitno,97) PDFmember,'LHAPDF set'
c
c      write(unitno,*)
c      write(unitno,*)
c     . '( [Jet definition and event cuts] )'
c      write(unitno,99) dsqrt(wsqmin),'m34min'
c      write(unitno,99) dsqrt(wsqmax),'m34max'
c      write(unitno,99) dsqrt(bbsqmin),'m56min'
c      write(unitno,99) dsqrt(bbsqmax),'m56max'
c      write(unitno,98) inclusive,'inclusive'
c      write(unitno,96) algorithm,'algorithm'
c      write(unitno,99) ptjetmin,'ptjetmin'
c      write(unitno,99) etajetmax,'etajetmax'
c      write(unitno,99) Rcut,'Rcut'
c      write(unitno,98) makecuts,'makecuts'
c      write(unitno,99) leptpt,'leptpt'
c      write(unitno,99) leptrap,'leptrap'
c      write(unitno,99) misspt,'misspt'
c      write(unitno,99) leptpt2,'leptpt2'
c      write(unitno,99) leptrap2,'leptrap2'
c      write(unitno,99) Rjlmin,'Rjlmin'
c      write(unitno,99) Rllmin,'Rllmin'
c      write(unitno,99) delyjjmin,'delyjjmin'
c      write(unitno,98) jetsopphem,'jetsopphem'
c      write(unitno,97) lbjscheme,'lbjscheme'
c      write(unitno,99) gammpt,'gammpt'
c      write(unitno,99) gammrap,'gammrap'
c      write(unitno,99) gammcone,'gammcone'
c      write(unitno,99) gammcut,'gammcut'
c
c      write(unitno,*)
c      write(unitno,*)
c     . '( [Anomalous couplings of the W and Z] )'
c      write(unitno,99) delg1_z,'delg1_z'
c      write(unitno,99) delk_z,'delk_z'
c      write(unitno,99) delk_g,'delk_g'
c      write(unitno,99) lambda_z,'lambda_z'
c      write(unitno,99) lambda_g,'lambda_g'
c      write(unitno,99) tevscale,'tevscale'
c
c      write(unitno,*)
c      write(unitno,*) 
c     . '( [How to resume/save a run] )'
c      write(unitno,98) readin,'readin'
c      write(unitno,98) writeout,'writeout'
c      write(unitno,96) ingridfile,'ingridfile'
c      write(unitno,96) outgridfile,'outgridfile'
c
c      write(unitno,99)

      return

c--- 55 format
   55 format(a20,f24.10,a4,f24.10,a1)
c--- 56 character format
   56 format('( PDF error set ',i3,'  --->',f13.3,' fb  )')
c--- 95 character format
   95 format(' (',5x,a9,' |',f18.5,f8.2,'%')
c--- 96 character format      
   96 format(' (',a20,12x,'[',a,']',' )')  
c--- 97 integer format      
   97 format(' (',i20,12x,'[',a,']',' )')  
c--- 98 logical format      
   98 format(' (',L20,12x,'[',a,']',' )')  
c--- 99 floating point format
   99 format(' (',f20.4,12x,'[',a,']',' )')  
      
      end
      
      
