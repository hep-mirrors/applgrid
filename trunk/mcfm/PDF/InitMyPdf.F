       subroutine InitMyPdf(name,iset)
        implicit none
        byte    name(*)     ! intent(in)

        integer iset
C       external         evolvePDF
        external         mypdfevolve
C        double precision x,Q2,pq(-6:6)
C        integer iloop,nf
        integer iflv,nloop
        double precision dy
C        double precision Q,ourpdf(-6:6),lhapdf(-6:6)
        integer nev
C        data nev/0/

        logical debug
        data debug/.false./

        integer maxlen 
        parameter(maxlen=200)
        character*(maxlen) f77name
        integer i,hitend

         write(6,*) '**Initpdf: starting set ',iset
              ! start the dglap evolution/convolution package
C         dy    = 0.01d0     ! 
          dy    = 0.1d0     ! the internal grid spacing (smaller->higher accuarcy)
                           ! 0.1 should provide at least 10^{-3} accuracy
         nloop = 2         ! the number of loops to initialise (max=3!)
         call hoppetStart(dy, nloop)

         hitend = 0
         do i = 1, maxlen
          if (name(i) .eq. 0) then
            hitend = i-1
            exit
          end if
          f77name(i:i) = achar(name(i))
 
         end do
         if (hitend == 0) then
          write(0,*) "Error in initPDFsetC: name is too long (max",
     $         maxlen," characters)"
          stop
         end if
         write(6,*) i,'InitMyPdf: name= ',f77name(1:hitend)


         ! initialise an LHAPDF set
C         write(6,*) ' init with: ',f77name(1:hitend)

         call InitPDFset(f77name(1:hitend))
C         call InitPDFsetByName(f77name(1:hitend))

         call InitPDF(iset)

         ! initialise our PDF using the LHAPDF subroutine for PDF-access
         ! (any other subroutine with same interface can be used in its place)
C         call hoppetAssign(evolvePDF)
         call hoppetAssign(mypdfevolve)

         ! for testing timing.
         !do i = 1, 100
         !   write(0,*) "hello"
         !   call hoppetAssign(mypdfevolve)
         !end do

       return 
       end   
