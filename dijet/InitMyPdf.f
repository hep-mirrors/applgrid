       subroutine InitMyPdf(name,iset)
         implicit none
         byte    name(*)        ! intent(in)

         integer iset

         integer maxlen 
         parameter(maxlen=200)
         character*(maxlen) f77name
         integer i,hitend
         
         write(6,*) '**Initpdf: starting set ',iset

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
     $           maxlen," characters)"
            stop
         end if

         write(6,*) i,'InitMyPdf: name= ',f77name(1:hitend)
         
C------  initialise an LHAPDF set
C        write(6,*) ' init with: ',f77name(1:hitend)
         
         call InitPDFset(f77name(1:hitend))
C        call InitPDFsetByName(f77name(1:hitend))

         call InitPDF(iset)

         return 
       end   
c
c----------------------------------------------------------------------
c      
       subroutine InitMyPdfSet(name)
         implicit none
         byte    name(*)        ! intent(in)

c         integer iset

         integer maxlen 
         parameter(maxlen=200)
         character*(maxlen) f77name
         integer i,hitend
         
c         write(6,*) '**Initpdf: starting set ',iset

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
     $           maxlen," characters)"
            stop
         end if

         write(6,*) i,'InitMyPdf: name= ',f77name(1:hitend)
         
C------  initialise an LHAPDF set
C        write(6,*) ' init with: ',f77name(1:hitend)
         
         call InitPDFset(f77name(1:hitend))
C        call InitPDFsetByName(f77name(1:hitend))

c         call InitPDF(iset)

         return 
       end   
c
c--------------------------------------------------------------------------------------
c
       subroutine InitMyPdfMember(iset)
         implicit none
         integer iset

         write(6,*) '**Initpdf: starting set ',iset

         call InitPDF(iset)

         return 
       end   


