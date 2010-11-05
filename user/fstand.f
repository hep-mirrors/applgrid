

C-----------------------------------------
C     dummy routines that should call 
C     the pdf and alphas routines
C-----------------------------------------
      double precision function fnalphas(Q)
      double precision Q      
      double precision fnalphas
      fnalphas = 1 
      return
      end

      subroutine fnpdf(x, Q, xf)
      double precision x, Q      
      double precision xf(13)
      do i=1, 13
         xf(i) = 1
      end do
      return
      end
C-----------------------------------------


      subroutine fstand
      implicit none

      double precision xsec(100)

      integer igrid
      integer idata

      integer Nbins 

      integer nloops

      integer igrids(100)

      integer inewgrids(100)

      integer n

      integer getnbins

      nloops = 1

C---  read a fastnlo grid -------------
C---  NB: pass in a vector of integers which returns the 
C---  ids of the grids which have been read in 
      call readfastnlogrids(igrids, "fnl0004.tab"//char(0))

C---  read native applgrid grids -------------
C---  NB: pass in an integer which returns the id of the grid
      call readgrid(igrids(1), "atlas-incljets04-eta1.root"//char(0))
      call readgrid(igrids(2), "atlas-incljets04-eta2.root"//char(0))
      call readgrid(igrids(3), "atlas-incljets04-eta3.root"//char(0))
      call readgrid(igrids(4), "atlas-incljets04-eta4.root"//char(0))
      call readgrid(igrids(5), "atlas-incljets04-eta5.root"//char(0))

C---  find out how many grids have been setup ----
      call ngrids(n) 
      call gridids( igrids ) 

C---  do the convolutions ---
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

