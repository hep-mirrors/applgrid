      subroutine nplotter_Wbbmas(p,wt,wt2,switch)
c--- Variable passed in to this routine:
c
c---      p:  4-momenta of particles in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)
c
c---     wt:  weight of this event
c
c---    wt2:  weight^2 of this event
c
c--- switch:  an integer equal to 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation
      implicit none
      include 'vegas_common.f'
      include 'constants.f'
      include 'histo.f'
      include 'jetlabel.f'
      integer ilept,inu,ibq,iba,ilt
      double precision p(mxpart,4),wt,wt2
      double precision yrap,pt,r,yraptwo,pttwo
      double precision ylept,ptlept,misset,yw,ptw,mw,ybq,ptbq,yba,ptba,
     & ylt,ptlt,yj1,ptj1,yj2,ptj2,yj3,ptj3,mbb,mjj,Rjl
      integer switch,n,nplotmax,nproc
      character*4 tag
      logical first,creatent,dswhisto
      common/outputflags/creatent,dswhisto
      common/nplotmax/nplotmax
      common/nproc/nproc
      data first/.true./
      save first
  
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

c--- Initialize dummy values for all quantities that could be plotted
      ylept=99d0
      yw=99d0
      ybq=99d0
      yba=99d0
      ylt=99d0
      yj1=99d0
      yj2=99d0
      yj3=99d0
      ptlept=-1d0
      misset=-1d0
      ptw=-1d0
      ptbq=-1d0
      ptba=-1d0
      ptlt=-1d0
      ptj1=-1d0
      ptj2=-1d0
      ptj3=-1d0
      mw=-1d0
      mbb=-1d0
      mjj=-1d0
      
      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag='book'
        goto 99
      else
c--- Add event in histograms
        tag='plot'
      endif

c--- Book and fill ntuple if that option is set, remembering to divide
c--- by # of iterations now that is handled at end for regular histograms
      if (creatent .eqv. .true.) then
        call bookfill(tag,p,wt/dfloat(itmx))  
      endif

************************************************************************
*                                                                      *
*     DEFINITIONS OF QUANTITIES TO PLOT                                *
*                                                                      *
*     Relevant processes are:                                          *
*                                                                      *
* 20  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +b(p5)+b~(p6) [massive]   *
* 25  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) +b(p5)+b~(p6) [massive]  *
*                                                                      *
************************************************************************

      if     ((nproc .eq. 20) .or. (nproc .eq. 21)
     &   .or. (nproc .eq. 420) .or. (nproc .eq. 421)
     &   .or. (nproc .eq. 422)) then
        ilept=4
	inu=3
      elseif ((nproc .eq. 25) .or. (nproc .eq. 26)
     &   .or. (nproc .eq. 425) .or. (nproc .eq. 426)
     &   .or. (nproc .eq. 427)) then
        ilept=3
	inu=4
      else
        write(6,*) 'Unanticipated process in nplotter_Wbbmas.f'
	stop 
      endif
      
      ylept=yrap(ilept,p)
      ptlept=pt(ilept,p)
      misset=pt(inu,p)
      
      yw=yraptwo(3,4,p)
      ptw=pttwo(3,4,p)
      mw=dsqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     &        -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)


c--- find the locations of the b,anti-b and light (or gluon) quarks      
      ibq=-1
      iba=-1
      ilt=-1
      if     (jetlabel(1) .eq. 'bq') then
        ibq=5
      elseif ((jetlabel(2) .eq. 'bq') .and. (jets .ge. 2)) then
        ibq=6
      elseif ((jetlabel(3) .eq. 'bq') .and. (jets .ge. 3)) then
        ibq=7
      endif
      
      if     (jetlabel(1) .eq. 'ba') then
        iba=5
      elseif ((jetlabel(2) .eq. 'ba') .and. (jets .ge. 2)) then
        iba=6
      elseif ((jetlabel(3) .eq. 'ba') .and. (jets .ge. 3)) then
        iba=7
      endif
      
      if     (jetlabel(1) .eq. 'pp') then
        ilt=5
      elseif ((jetlabel(2) .eq. 'pp') .and. (jets .ge. 2)) then
        ilt=6
      elseif ((jetlabel(3) .eq. 'pp') .and. (jets .ge. 3)) then
        ilt=7
      endif
      
c--- this case only occurs for processes 421, 426
      if     (jetlabel(1) .eq. 'bb') then
        ibq=5
      elseif ((jetlabel(2) .eq. 'bb') .and. (jets .ge. 2)) then
        ibq=6
      endif
      
      Rjl=99d0
      if (ibq .gt. 0) then
        ybq=yrap(ibq,p)
        ptbq=pt(ibq,p)
	Rjl=min(Rjl,R(p,ilept,ibq))
      endif
      if (iba .gt. 0) then
        yba=yrap(iba,p)
        ptba=pt(iba,p)
	Rjl=min(Rjl,R(p,ilept,iba))
      endif
      if (ilt .gt. 0) then
        ylt=yrap(ilt,p)
        ptlt=pt(ilt,p)
	Rjl=min(Rjl,R(p,ilept,ilt))
      endif
       
      yj1=yrap(5,p)
      ptj1=pt(5,p)
      if (jets .ge. 2) then
        yj2=yrap(6,p)
        ptj2=pt(6,p)
      endif
      if (jets .ge. 3) then
        yj3=yrap(6,p)
        ptj3=pt(6,p)
      endif

      if (jets .gt. 1) then
        mjj=dsqrt((p(5,4)+p(6,4))**2-(p(5,1)+p(6,1))**2
     &           -(p(5,2)+p(6,2))**2-(p(5,3)+p(6,3))**2)
      endif
      if ((ibq .gt. 0) .and. (iba .gt. 0)) then
        mbb=dsqrt((p(ibq,4)+p(iba,4))**2-(p(ibq,1)+p(iba,1))**2
     &           -(p(ibq,2)+p(iba,2))**2-(p(ibq,3)+p(iba,3))**2)
      endif
      
c      write(6,*) 'In Wbbmas, switch=',switch
c      write(6,*) 'In Wbbmas, jets=',jets
c      write(6,*) 'In Wbbmas, jetlabel=',jetlabel
c      write(6,*)
      
************************************************************************
*                                                                      *
*     FILL HISTOGRAMS                                                  *
*                                                                      *
************************************************************************

c--- Call histogram routines
   99 continue

c--- "n" will count the number of histograms
      n=1              

c--- Syntax of "bookplot" routine is:
c
c---   call bookplot(n,tag,titlex,var,wt,wt2,xmin,xmax,dx,llplot)
c
c---        n:  internal number of histogram
c---      tag:  "book" to initialize histogram, "plot" to fill
c---   titlex:  title of histogram
c---      var:  value of quantity being plotted
c---       wt:  weight of this event (passed in)
c---      wt2:  weight of this event (passed in)
c---     xmin:  lowest value to bin
c---     xmax:  highest value to bin
c---       dx:  bin width
c---   llplot:  equal to "lin"/"log" for linear/log scale

      call bookplot(n,tag,'Lepton rapidity',ylept,wt,wt2,
     &              -5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'Lepton pt',ptlept,wt,wt2,
     &              0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'Missing Et',misset,wt,wt2,
     &              0d0,200d0,5d0,'log')
      n=n+1

      call bookplot(n,tag,'W boson rapidity',yw,wt,wt2,
     &              -5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'W boson pt',ptw,wt,wt2,
     &              0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'W boson invariant mass',mw,wt,wt2,
     &              0d0,200d0,5d0,'lin')
      n=n+1

      call bookplot(n,tag,'b-quark jet rapidity',ybq,wt,wt2,
     &              -5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'b-quark jet pt',ptbq,wt,wt2,
     &              0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'anti-b-quark jet rapidity',yba,wt,wt2,
     &              -5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'anti-b-quark jet pt',ptba,wt,wt2,
     &              0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'(b,anti-b) jet invariant mass',mbb,wt,wt2,
     &              0d0,800d0,10d0,'log')
      n=n+1

      call bookplot(n,tag,'Jet 1 rapidity',yj1,wt,wt2,
     &              -5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 1 pt',ptj1,wt,wt2,
     &              0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'Jet 2 rapidity',yj2,wt,wt2,
     &              -5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 2 pt',ptj2,wt,wt2,
     &              0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'(jet 1,jet 2) invariant mass',mjj,wt,wt2,
     &              0d0,800d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'Minimum R(lepton,jet)',Rjl,wt,wt2,
     &              0d0,6d0,0.1d0,'lin')
      n=n+1

      call bookplot(n,tag,'light jet rapidity',ylt,wt,wt2,
     &              -5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'light jet pt',ptlt,wt,wt2,
     &              0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'Jet 3 rapidity',yj3,wt,wt2,
     &              -5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 3 pt',ptj3,wt,wt2,
     &              0d0,200d0,5d0,'log')
      n=n+1

c      call bookplot(n,tag,'DeltaRe5',re5,wt,wt2,0d0,5d0,0.1d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'y5',y5,wt,wt2,-yjet,yjet,0.2d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'pt5',pt5,wt,wt2,0d0,ptjet,2d0,'lin')
c      n=n+1

  
************************************************************************
*                                                                      *
*     FINAL BOOKKEEPING                                                *
*                                                                      *
************************************************************************

c--- We have over-counted the number of histograms by 1 at this point
      n=n-1

c--- Ensure the built-in maximum number of histograms is not exceeded    
      call checkmaxhisto(n)

c--- Set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      
      return 
      end
      
