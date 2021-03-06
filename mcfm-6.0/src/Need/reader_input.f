      subroutine reader_input
************************************************************************
*   Routine to read in the file input.DAT, which is a consolidated     *
*   form of all the input files and new to version 3.4 of MCFM         *
************************************************************************
      implicit none
      include 'constants.f'
      include 'debug.f'
      include 'new_pspace.f'
      include 'virtonly.f'
      include 'realonly.f'
      include 'noglue.f'
      include 'realwt.f'
      include 'lc.f'
      include 'cutoff.f'
      include 'maxwt.f'
      include 'masses.f'
      include 'process.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'zerowidth.f'
      include 'removebr.f'
      include 'flags.f'
      include 'clustering.f'
      include 'anomcoup.f'
      include 'gridinfo.f'
      include 'verbose.f'
      include 'limits.f'
      include 'workdir.f'
      include 'jetcuts.f'
      include 'leptcuts.f'
      include 'lhapdf.f'
      include 'alfacut.f'
      include 'betacut.f'
      include 'pdlabel.f'
      include 'qcdcouple.f'
      include 'nlooprun.f'
      include 'initialscales.f'
      include 'stopscales.f'
      include 'vanillafiles.f'
      include 'frag.f'

c     P.S. flag using grid
      include 'ptilde.f'
      include 'APPLinclude.f'
c
      character*72 inputfile,getinput
      character*90 line
      character*4 part
      character*30 runstring
      integer nargs,lenocc,lenarg
      logical spira
      logical creatent,dswhisto,dryrun,makecuts
      integer nmin,nmax,n2,n3
      integer nproc,ih1,ih2,itmx1,itmx2,ncall1,ncall2,idum,origij
      integer NPTYPE,NGROUP,NSET
      double precision rtsmin,sqrts,factor
      double precision mbbmin,mbbmax,Mwmin,Mwmax
      double precision Rcut
      logical technicalincluded
      double precision ran2,randummy
      double precision cmass,bmass
      double precision mass2,width2,mass3,width3
      double precision amz,alphas
      
C---- prototypes for argument passing -----
      integer cargc
      integer clenargv

      common/couple/amz
      
      common/breit/n2,n3,mass2,width2,mass3,width3
      
      common/spira/spira
      common/nmin/nmin
      common/nmax/nmax
      common/rtsmin/rtsmin
 
      common/outputflags/creatent,dswhisto      

      common/nproc/nproc
      common/part/part
      common/runstring/runstring
      common/energy/sqrts
      common/density/ih1,ih2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/ranno/idum
      common/dryrun/dryrun
      
      common/pdflib/NPTYPE,NGROUP,NSET
      
      common/Rcut/Rcut
      common/makecuts/makecuts

      common/qmass/cmass,bmass

      common/origij/origij

      save /ranno/
      
      verbose=.true.

c--- work out the name of the input file and open it
c--- use functions defined in main.cxx since this is 
c--- c++ linked now 

C     nargs=iargc()
C      nargs=iargc()
C      if (nargs .ge. 1) then
C        call getarg(1,inputfile)
C      else
C        inputfile='input.DAT'
C      endif
      
C     lenarg=lenocc(inputfile)

      nargs=cargc()

      write (6,*) "argc=", nargs

      lenarg = 0

      if (nargs .ge. 1) then
         inputfile='duffer'
         call cargv(1,inputfile)
         lenarg = clenargv(1)
C        call get_command_argument(1, inputfile, cmdlength, cmdstatus)
C        call getarg(1,inputfile)
         write (6,*) "inputfile=", inputfile, "lenarg=", lenarg
      else
        inputfile='Winput.DAT'
      endif
     

      if ((lenarg.lt.4).or.(inputfile(lenarg-3:lenarg).ne.'.DAT')) then
        workdir=inputfile
c--- truncate if the directory / is included
        if (workdir(lenarg:lenarg) .eq. '/') then
           workdir(lenarg:lenarg)=' '
           lenarg=lenarg-1
        endif
        if (nargs .ge. 2) then
          call getarg(2,getinput)
          inputfile=workdir(1:lenarg)//'/'//getinput
        else  
          inputfile=workdir(1:lenarg)//'/input.DAT'
        endif
      else
        workdir=''
      endif
            
      write(6,*) '* Using input file named ',inputfile

      open(unit=20,file=inputfile,status='old',err=999)
      call checkversion(20,inputfile)

c--- read-in the user inputs

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- flags for the mode of MCFM   
      read(20,*) evtgen
      if (verbose) call writeinput(6,' * ',' ','evtgen')
      read(20,*) creatent
      if (verbose) call writeinput(6,' * ',' ','creatent')
      read(20,*) skipnt
      if (verbose) call writeinput(6,' * ',' ','skipnt')
      read(20,*) dswhisto
      if (verbose) call writeinput(6,' * ',' ','dswhisto')
c
c     P.S.  added to read flag for grids
c
       read(20,*) creategrid
       if (verbose) call writeinput(6,' * ',' ','creategrid')
cc
cc
      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- general options 
      read(20,*) nproc
      if (verbose) call writeinput(6,' * ',' ','nproc')
      read(20,*) part
      if (verbose) call writeinput(6,' * ',' ','part')
      read(20,*) runstring
      if (verbose) call writeinput(6,' * ',' ','runstring')
      read(20,*) sqrts
      if (verbose) call writeinput(6,' * ',' ','sqrts')
      read(20,*) ih1
      if (verbose) call writeinput(6,' * ',' ','ih1')
      read(20,*) ih2
      if (verbose) call writeinput(6,' * ',' ','ih2')
      read(20,*) hmass
      if (verbose) call writeinput(6,' * ',' ','hmass')
      read(20,*) scale
      if (verbose) call writeinput(6,' * ',' ','scale')
      read(20,*) facscale
      if (verbose) call writeinput(6,' * ',' ','facscale')
      
      initrenscale_L=0d0
      initfacscale_L=0d0
      initrenscale_H=0d0
      initfacscale_H=0d0
c--- catch special scale choices for stop+b process
      if (((nproc .eq. 231) .or. (nproc .eq. 236)      
     . .or.(nproc .eq. 241) .or. (nproc .eq. 246)      
     . .or.(nproc .eq. 242) .or. (nproc .eq. 247)) .and.      
     .    (scale .eq. 0d0) .and. (facscale .eq. 0d0)) then
        read(20,*) initrenscale_L
	renscale_L=initrenscale_L
        if (verbose) call writeinput(6,' * ',' ','renscale_L')
        read(20,*) initfacscale_L
	facscale_L=initfacscale_L
        if (verbose) call writeinput(6,' * ',' ','facscale_L')
        read(20,*) initrenscale_H
	renscale_H=initrenscale_H
        if (verbose) call writeinput(6,' * ',' ','renscale_H')
        read(20,*) initfacscale_H
	facscale_H=initfacscale_H
        if (verbose) call writeinput(6,' * ',' ','facscale_H')
	scale=initrenscale_H
	facscale=initfacscale_H        
      endif

      read(20,*) dynamicscale
      if (verbose) call writeinput(6,' * ',' ','dynamicscale')
      read(20,*) zerowidth
      if (verbose) call writeinput(6,' * ',' ','zerowidth')
      read(20,*) removebr
      if (verbose) call writeinput(6,' * ',' ','removebr')
      read(20,*) itmx1
      if (verbose) call writeinput(6,' * ',' ','itmx1')
      read(20,*) ncall1
      if (verbose) call writeinput(6,' * ',' ','ncall1')
      read(20,*) itmx2
      if (verbose) call writeinput(6,' * ',' ','itmx2')
      read(20,*) ncall2
      if (verbose) call writeinput(6,' * ',' ','ncall2')
      read(20,*) origij
      if (verbose) call writeinput(6,' * ',' ','ij')
      read(20,*) dryrun
      if (verbose) call writeinput(6,' * ',' ','dryrun')
      read(20,*) Qflag
      if (verbose) call writeinput(6,' * ',' ','Qflag')
      read(20,*) Gflag
      if (verbose) call writeinput(6,' * ',' ','Gflag')
      
      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- heavy quark masses 
      read(20,*) mt
      if (verbose) call writeinput(6,' * ',' ','top mass')
      read(20,*) mb
      if (verbose) call writeinput(6,' * ',' ','bottom mass')
      read(20,*) mc
      if (verbose) call writeinput(6,' * ',' ','charm mass')

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- pdf options 
      read(20,*) pdlabel
      if (verbose) call writeinput(6,' * ',' ','pdlabel')
      read(20,*) NGROUP
      if (verbose) call writeinput(6,' * ',' ','NGROUP')
      read(20,*) NSET
      if (verbose) call writeinput(6,' * ',' ','NSET')
      read(20,*) PDFname
      if (verbose) call writeinput(6,' * ',' ','LHAPDF group')
      read(20,*) PDFmember
      if (verbose) call writeinput(6,' * ',' ','LHAPDF set')

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- jets and cuts options 
      read(20,*) Mwmin
      if (verbose) call writeinput(6,' * ',' ','m34min')
      read(20,*) Mwmax 
      if (verbose) call writeinput(6,' * ',' ','m34max')
      read(20,*) mbbmin
      if (verbose) call writeinput(6,' * ',' ','m56min')
      read(20,*) mbbmax 
      if (verbose) call writeinput(6,' * ',' ','m56max')
      read(20,*) inclusive
      if (verbose) call writeinput(6,' * ',' ','inclusive')
      read(20,*) algorithm
      if (verbose) call writeinput(6,' * ',' ','algorithm')
      read(20,*) ptjetmin
      if (verbose) call writeinput(6,' * ',' ','ptjetmin')
      read(20,*) etajetmin 
      if (verbose) call writeinput(6,' * ',' ','etajetmin')
      read(20,*) etajetmax 
      if (verbose) call writeinput(6,' * ',' ','etajetmax')
      read(20,*) Rcut
      if (verbose) call writeinput(6,' * ',' ','Rcut')
      read(20,*) makecuts
      if (verbose) call writeinput(6,' * ',' ','makecuts')
      read(20,*) leptpt
      if (verbose) call writeinput(6,' * ',' ','leptpt')
      read(20,*) leptrap
      if (verbose) call writeinput(6,' * ',' ','leptrap')
      read(20,*) misspt
      if (verbose) call writeinput(6,' * ',' ','misspt')
      read(20,*) leptpt2
      if (verbose) call writeinput(6,' * ',' ','leptpt2')
      read(20,*) leptrap2
      if (verbose) call writeinput(6,' * ',' ','leptrap2')
      read(20,*) mtrans34cut
      if (verbose) call writeinput(6,' * ',' ','mtrans34cut')
      read(20,*) Rjlmin
      if (verbose) call writeinput(6,' * ',' ','Rjlmin')
      read(20,*) Rllmin
      if (verbose) call writeinput(6,' * ',' ','Rllmin')
      read(20,*) delyjjmin
      if (verbose) call writeinput(6,' * ',' ','delyjjmin')
      read(20,*) jetsopphem 
      if (verbose) call writeinput(6,' * ',' ','jetsopphem')
      read(20,*) lbjscheme 
      if (verbose) call writeinput(6,' * ',' ','lbjscheme')
      read(20,*) ptbjetmin
      if (verbose) call writeinput(6,' * ',' ','ptbjetmin')
      read(20,*) etabjetmax
      if (verbose) call writeinput(6,' * ',' ','etabjetmax')

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- settings for photon processes 
      read(20,*) frag
      if (verbose) call writeinput(6,' * ',' ','frag')
      read(20,*) fragset
      if (verbose) call writeinput(6,' * ',' ','fragset')
      read(20,*) frag_scale
      if (verbose) call writeinput(6,' * ',' ','frag_scale')
      read(20,*) gammpt
      if (verbose) call writeinput(6,' * ',' ','gammpt')
      read(20,*) gammrap
      if (verbose) call writeinput(6,' * ',' ','gammrap')
      read(20,*) Rgalmin 
      if (verbose) call writeinput(6,' * ',' ','Rgalmin')
      read(20,*) cone_ang
      if (verbose) call writeinput(6,' * ',' ','cone_ang')
      read(20,*) epsilon_h
      if (verbose) call writeinput(6,' * ',' ','epsilon_h')

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- anomalous couplings 
      read(20,*) delg1_z
      if (verbose) call writeinput(6,' * ',' ','delg1_z')
      read(20,*) delk_z
      if (verbose) call writeinput(6,' * ',' ','delk_z')
      read(20,*) delk_g
      if (verbose) call writeinput(6,' * ',' ','delk_g')
      read(20,*) lambda_z
      if (verbose) call writeinput(6,' * ',' ','lambda_z')
      read(20,*) lambda_g
      if (verbose) call writeinput(6,' * ',' ','lambda_g')
      read(20,*) tevscale
      if (verbose) call writeinput(6,' * ',' ','tevscale')

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- grid information 
      read(20,*) readin
      if (verbose) call writeinput(6,' * ',' ','readin')
      read(20,*) writeout
      if (verbose) call writeinput(6,' * ',' ','writeout')
      read(20,*) ingridfile
      if (verbose) call writeinput(6,' * ',' ','ingridfile')
      read(20,*) outgridfile
      if (verbose) call writeinput(6,' * ',' ','outgridfile')

      if (verbose) write(6,*)

c--- check if the contents of technical.DAT are included here
c--- (the default behaviour going forward)
      technicalincluded=.false.

      read(20,99,end=88) line
      read(20,99,end=88) line
c      write(6,*) line
      if (line(2:10) .ne. 'Technical') goto 88
      technicalincluded=.true.

   88 continue
      
c--- if the contents of technical.DAT are not included, close the input
c--- file and open technical.DAT instead; otherwise continue on
      if (technicalincluded .eqv. .false.) then
        close(20)
        open(unit=20,file='technical.DAT',status='old',err=999)
        call checkversion(20,'technical.DAT')
      endif

      if (verbose) write(6,*) '* [Technical parameters that'//
     .                        ' should not normally be changed]'
      if (verbose) write(6,*)

c---- read-in the technical parameters

      read(20,*) debug
      if (verbose) call writeinput(6,' * ',' ','debug')
      read(20,*) verbose
      if (verbose) call writeinput(6,' * ',' ','verbose')
      read(20,*) new_pspace
      if (verbose) call writeinput(6,' * ',' ','new_pspace')
      read(20,*) virtonly
      if (verbose) call writeinput(6,' * ',' ','virtonly')
      read(20,*) realonly
      if (verbose) call writeinput(6,' * ',' ','realonly')
      read(20,*) spira
      if (verbose) call writeinput(6,' * ',' ','spira')
      read(20,*) noglue
      if (verbose) call writeinput(6,' * ',' ','noglue')
      read(20,*) ggonly
      if (verbose) call writeinput(6,' * ',' ','ggonly')
      read(20,*) gqonly
      if (verbose) call writeinput(6,' * ',' ','gqonly')
      read(20,*) vanillafiles
      if (verbose) call writeinput(6,' * ',' ','vanillafiles')
      read(20,*) nmin
      if (verbose) call writeinput(6,' * ',' ','nmin')
      read(20,*) nmax
      if (verbose) call writeinput(6,' * ',' ','nmax')
      read(20,*) clustering
      if (verbose) call writeinput(6,' * ',' ','clustering')
      read(20,*) realwt
      if (verbose) call writeinput(6,' * ',' ','realwt')
      read(20,*) colourchoice
      if (verbose) call writeinput(6,' * ',' ','colourchoice')
      read(20,*) rtsmin
      if (verbose) call writeinput(6,' * ',' ','rtsmin')
      read(20,*) cutoff
      if (verbose) call writeinput(6,' * ',' ','cutoff')
      read(20,*) aii
      if (verbose) call writeinput(6,' * ',' ','aii')
      read(20,*) aif
      if (verbose) call writeinput(6,' * ',' ','aif')
      read(20,*) afi
      if (verbose) call writeinput(6,' * ',' ','afi')
      read(20,*) aff
      if (verbose) call writeinput(6,' * ',' ','aff')
      read(20,*) bfi
      if (verbose) call writeinput(6,' * ',' ','bfi')
      read(20,*) bff
      if (verbose) call writeinput(6,' * ',' ','bff')
      if (verbose) write(6,*)
      close(unit=20)

      if ((etajetmin .lt. 0d0) .or. (etajetmax .lt. 0d0)) then
        write(6,*) 'etajetmin and etajetmax are absolute values,'
	write(6,*) ' please reset to a positive value.'
	stop
      endif

c--- reset values of the alpha parameters, for specific runstrings
      if (runstring(1:5) .eq. 'alpha') then
        if     (runstring(6:9) .eq. '1111') then
          aii=1d0
          aif=1d0
          afi=1d0
          aff=1d0
        elseif (runstring(6:9) .eq. '0111') then
          aii=0.1d0
          aif=1d0
          afi=1d0
          aff=1d0
        elseif (runstring(6:9) .eq. '1011') then
          aii=1d0
          aif=0.1d0
          afi=1d0
          aff=1d0
        elseif (runstring(6:9) .eq. '1101') then
          aii=1d0
          aif=1d0
          afi=0.1d0
          aff=1d0
        elseif (runstring(6:9) .eq. '1110') then
          aii=1d0
          aif=1d0
          afi=1d0
          aff=0.1d0
        elseif (runstring(6:9) .eq. '0001') then
          aii=0.1d0
          aif=0.1d0
          afi=0.1d0
          aff=1d0
        elseif (runstring(6:9) .eq. '0000') then
          aii=0.1d0
          aif=0.1d0
          afi=0.1d0
          aff=0.1d0
        else
          write(6,*) 'runstring=alpha.... reserved for internal tests'
          stop
        endif
        write(6,*) 'Values of alpha parameters reset, for testing:'
        write(6,*) '  alpha_ii=',aii
        write(6,*) '  alpha_if=',aif
        write(6,*) '  alpha_fi=',afi
        write(6,*) '  alpha_ff=',aff
        write(6,*)
      endif

      if     (index(runstring,'noglue') .gt. 0) then
        noglue=.true.
      elseif (index(runstring,'ggonly') .gt. 0) then
        ggonly=.true.
      elseif (index(runstring,'gqonly') .gt. 0) then
        gqonly=.true.
      endif
      
      if     (index(runstring,'mc1.3') .gt. 0) then
        mc=1.3d0
	mcsq=mc**2
      elseif (index(runstring,'mc1.4') .gt. 0) then
        mc=1.4d0
	mcsq=mc**2
      elseif (index(runstring,'mc1.5') .gt. 0) then
        mc=1.5d0
	mcsq=mc**2
      endif
      
      if (runstring(1:3) .eq. 'mlm') then
        write(6,*) 'WARNING: cross sections divided by Ecm**2'
	write(6,*)
      endif
      
      if (noglue) then
        write(6,*) 'WARNING: no gluon contribution included in PDF'
	write(6,*)
      endif
      if (ggonly) then
        write(6,*) 'WARNING: only gluon-gluon flux included'
	write(6,*)
      endif
      if (gqonly) then
        write(6,*) 'WARNING: only gluon-quark flux included'
	write(6,*)
      endif
      
c-----initialize various quantities

c--- save initial scale choices (that may be changed later)
      initscale=scale
      initfacscale=facscale

c--- assign squared masses for b- and c-quarks
      if (abs(mb) .gt. 1d-8) then
        mbsq=mb**2
      else
        mbsq=4.75d0**2
      endif
      if (abs(mc) .gt. 1d-8) then
        mcsq=mc**2
      else
        mcsq=1.5d0**2
      endif
      
c--- set-up mass window cuts
      bbsqmin=mbbmin**2
      bbsqmax=mbbmax**2

      wsqmin=Mwmin**2
      wsqmax=Mwmax**2

c--- set-up the variables for the process we wish to consider
      call chooser

c--- set-up the random number generator with a negative seed
      idum=-abs(origij)
      randummy=ran2()

c--- initialize masses for alpha_s routine
      cmass=dsqrt(mcsq)
      bmass=dsqrt(mbsq)

c--- E-M gauge invariance requires that delg1_g=0
      delg1_g=0d0

c--- check that we have a valid value of 'part'
      if ( (part .ne. 'lord') .and. (part .ne. 'real') .and.
     .     (part .ne. 'virt') .and. (part .ne. 'tota') ) then
        if    ( (part .eq. 'todk') .and.
     .          ((case .eq. 'bq_tpq') .or. (case .eq. 't_bbar')
     .      .or. (case .eq. 'W_twdk')) ) then
c--- this is an allowed combination
        elseif ( (part .eq. 'frag') .and.
     .       ((case .eq. 'Wgamma') .or. (case .eq. 'Zgamma')
     &       .or.(case.eq.'gamgam'))) then
c--- this is an allowed combination
        else 
          write(6,*) 'part=',part,' is not a valid option'
          write(6,*) 'for this process number.'
          stop     
        endif
      endif      

c--- set up the default choices of static scale, if required
      if (scale .lt. 0d0) then
	if     (scale .eq. -2d0) then
	  factor=0.25d0
	elseif (scale .eq. -3d0) then
	  factor=0.5d0
	elseif (scale .eq. -4d0) then
	  factor=0.75d0
	elseif (scale .eq. -5d0) then
	  factor=1d0
	elseif (scale .eq. -6d0) then
	  factor=2d0
	elseif (scale .eq. -7d0) then
	  factor=4d0
        else
	  factor=1d0
	endif	  
        if ((n2+n3 .ne. 0) .or. (case .eq. 'tt_tot')) then
        scale=factor*(dfloat(n2)*mass2+dfloat(n3)*mass3)/dfloat(n2+n3)
c--- special cases where Higgs mass is neither mass2 nor mass3
        if ((case(1:1) .eq. 'H') .or. (case .eq. 'WHbbar')
     . .or. (case .eq. 'ZHbbar') .or. (case .eq. 'qq_Hqq')) then
          scale=factor*hmass
	endif
c--- special case for t-tbar production
        if (case .eq. 'tt_tot') scale=factor*mt	
        as=alphas(scale,amz,nlooprun)
        ason2pi=as/twopi
        ason4pi=as/fourpi
        gsq=fourpi*as
        musq=scale**2
        write(6,*)
        write(6,*)'************* Strong coupling, alpha_s  ************'
        write(6,*)'*                                                  *'
        write(6,49)'alpha_s (scale)',gsq/fourpi
        write(6,49)'alpha_s (zmass)',amz
        write(6,50)' (using ',nlooprun,'-loop running of alpha_s)'  
        write(6,*)'****************************************************'
        write(6,*)
        write(6,*)'****************************************************'
        write(6,76) scale
        write(6,*)'****************************************************'
        else
        write(6,*) 'Invalid choice of renormalization scale!'
        stop
        endif
      endif
      if (facscale .lt. 0d0) then
	if     (facscale .eq. -2d0) then
	  factor=0.25d0
	elseif (facscale .eq. -3d0) then
	  factor=0.5d0
	elseif (facscale .eq. -4d0) then
	  factor=0.75d0
	elseif (facscale .eq. -5d0) then
	  factor=1d0
	elseif (facscale .eq. -6d0) then
	  factor=2d0
	elseif (facscale .eq. -7d0) then
	  factor=4d0
        else
	  factor=1d0
	endif	  
        if ((n2+n3 .ne. 0) .or. (case .eq. 'tt_tot')) then
        facscale=factor*
     .           (dfloat(n2)*mass2+dfloat(n3)*mass3)/dfloat(n2+n3)
c--- special cases where Higgs mass is neither mass2 nor mass3
        if ((case(1:1) .eq. 'H') .or. (case .eq. 'WHbbar')
     . .or. (case .eq. 'ZHbbar') .or. (case .eq. 'qq_Hqq')) then
          facscale=factor*hmass
	endif
c--- special case for t-tbar production
        if (case .eq. 'tt_tot') facscale=factor*mt	
        write(6,*)
        write(6,*)'****************************************************'
        write(6,77) facscale
        write(6,*)'****************************************************'
       else
        write(6,*) 'Invalid choice of factorization scale!'
        stop
        endif
      endif
      
      return

   49 format(' *  ',a20,f12.8,16x,'*')
   50 format(' *  ',6x,a8,i1,a25,8x,'*')
   76 format(' *      Renormalization scale =',f7.2,'              *')   
   77 format(' *        Factorization scale =',f7.2,'              *')   
   99 format(a90)

  999 continue
      write(6,*) 'Problem reading ',inputfile
      write(6,*)
      write(6,*) 'Refer to documentation for the format of input.DAT'
      write(6,*)
      stop

      end
      
