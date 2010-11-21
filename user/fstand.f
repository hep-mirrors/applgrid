

C----------------------------------------------------------
C     dummy routines that call  the pdf and alphas routines
C     this example, they just call the lhapdf routine
C     for apdf fit, call your own routines 
C----------------------------------------------------------
      double precision function fnalphas(Q)
      double precision Q      
      double precision alphaspdf
      fnalphas = alphaspdf(Q)
      return
      end

      subroutine fnpdf(x, Q, xf)
      double precision x, Q      
      double precision xf(13)
      call evolvePDF(x, Q, xf)
      return
      end
C----------------------------------------------------------


      subroutine fstand
      implicit none

      double precision xsec(100)

      integer idata

C---  counter of grid ids - this is fortran, so need ---
C---  to differentiate them somehow                  ---      
      integer igrid
      integer igrids(100)
      
      
      integer Nbins 
      integer inewgrids(100)
      integer n
      integer getnbins

C---  lhapdf set       
      integer iset

      integer dummy

C---  first set up the your pdfs however you like -----

      iset = 0

      call initPDFSetByName("cteq6mE.LHgrid") 
      call initpdf(iset)

C--- next simply read in the grids ----

C---  read a fastnlo grid -------------
C---  NB: pass in a vector of integers to be be filled with the 
C---  ids of the grids which have been read in 
C---  also note the need for the terminating //char(0) for 
C---  strings passed to C code
      call readfastnlogrids( igrids, "fnl0004.tab"//char(0) )

      call ngrids(n) 
      
      write(6,*) n, " grids in ", "fnl0004.tab" 

C---  read native applgrid grids -------------

C---  NB: pass in an integer which will be filled with the id 
C---  of the grid read in

C---  Note, we ignore the values returned in dummy here,
C---  since we will get back all the identifies in the call 
C---  to gridids() later. This would not be the case if 
C---  we wanted to actually know which grid was for which
C---  cross section  
      call readgrid( dummy, "atlas-incljets04-eta1.root"//char(0) )
      call readgrid( dummy, "atlas-incljets04-eta2.root"//char(0) )
      call readgrid( dummy, "atlas-incljets04-eta3.root"//char(0) )
      call readgrid( dummy, "atlas-incljets04-eta4.root"//char(0) )
      call readgrid( dummy, "atlas-incljets04-eta5.root"//char(0) )

C---  find out how many grids have been setup ----
      call ngrids(n) 
C---  and get all their identiers ---
      call gridids( igrids ) 

      write(6,*) n, "grids in total"

C---  now do the convolutions ---

      do igrid=1, n

C---     how many bins in this cross section ---
         Nbins = getnbins(igrids(igrid))

         write(6,*) "grid id ", igrids(igrid), " with ", Nbins, "bins"

C---     print out the documentation for each grid --- 
C        call printgriddoc(igrids(igrid))

C---     convolute this grid ---
         call convolute(igrids(igrid), xsec)
         
C---     print out the results ---
         do idata=1, Nbins
            write(6,*) "xsec(", idata, ")=", xsec(idata)
         end do

      end do
      
C---  free the grid storage ---
      call releasegrids

      end

