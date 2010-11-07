

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
C-----------------------------------------


      subroutine fstand
      implicit none

      double precision xsec(100)

      integer idata

C---  counter of grid ids - this is fortran, so need ---
C---  to differentiate them somehow                  ---      
      integer igrid
      integer igrids(100)
      
      
      integer Nbins 
      integer nloops
      integer inewgrids(100)
      integer n
      integer getnbins

C---  lhapdf set       
      integer iset

      nloops = 1

C---  first set up the your pdfs however you like -----

      iset = 0

      call initPDFSet("/usr/local/share/lhapdf/PDFsets/cteq6mE.LHgrid") 
      call initpdf(iset)

C--- next simply read in the grids ----

C---  read a fastnlo grid -------------
C---  NB: pass in a vector of integers which returns the 
C---  ids of the grids which have been read in 
      call readfastnlogrids(igrids, "fnl0004.tab"//char(0))

      call ngrids(n) 
      
      write(6,*) n, " grids in ", "fnl0004.tab" 

C---  read native applgrid grids -------------
C---  NB: pass in an integer which returns the id of the grid
      call readgrid(igrid, "atlas-incljets04-eta1.root"//char(0))

      call readgrid(igrid, "atlas-incljets04-eta2.root"//char(0))
      call readgrid(igrid, "atlas-incljets04-eta3.root"//char(0))
      call readgrid(igrid, "atlas-incljets04-eta4.root"//char(0))
      call readgrid(igrid, "atlas-incljets04-eta5.root"//char(0))

C---  find out how many grids have been setup ----
      call ngrids(n) 
      call gridids( igrids ) 

      write(6,*) n, "grids in total"


C---  now do the convolutions ---

      do igrid=1, n

C---     how many bins in this cross section ---
         Nbins = getnbins(igrids(igrid))

         write(6,*) "grid ", igrids(igrid), Nbins, nloops

C---     convolute this crid ---
C---     NB: here we have dummy, nonsensical pdf and alphas
C---         routines, so the actual values are not important 
         call convolute(igrids(igrid), nloops, xsec)
         
C---     print out the results ---
         do idata=1, Nbins
            write(6,*) "xsec(", idata, ")=", xsec(idata)
         end do

      end do
      
C---  free the grid storage ---
      call releasegrids

      end

