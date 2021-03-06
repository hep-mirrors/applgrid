*****************
* LHAPDF version*
*****************
      subroutine pdfwrap
      implicit none
      include 'masses.f'
      include 'lhapdf.f'
      include 'nlooprun.f'
      include 'PDFerrors.f'
      include 'pdlabel.f'
      double precision amz,alphasPDF
      logical validPDF
      character*30 oldPDFname
      integer i,iorder

      common/couple/amz

      
c      if (newinput .eqv. .false.) then
c        open(unit=21,file='lhapdf.DAT',status='old',err=999)
c        call checkversion(21,'lhapdf.DAT')
c        read(21,*) PDFname
c        read(21,*) PDFmember            
c        close(21)
c      endif
      
      oldPDFname=PDFname
      validPDF=.false.
      i=0
   20 continue
      i=i+1    
      if ((oldPDFname(i:i) .eq. '.') .or.
     .    (oldPDFname(i:i) .eq. ' ') .or.
     .    (oldPDFname(i:i) .eq. '[')) then
        validPDF=.true.
        if (oldPDFname(i:i+6) .eq. '.LHgrid') then        
          PDFname=oldPDFname(1:i-1)//'.LHgrid'
        else
          PDFname=oldPDFname(1:i-1)//'.LHpdf'
        endif
      endif  
      if ((i .lt. 20) .and. (validPDF .eqv. .false.)) goto 20
      
      if (validPDF .eqv. .false.) then
        write(6,*) 'Problem with PDFname'
        write(6,*)
        stop
      endif
      
      write(6,*)
      write(6,*) '*******************************************'
      write(6,*) '*     MCFM is calling LHAPDF              *'
      write(6,*) '*                                         *'
      write(6,98) 'PDFname',PDFname(1:29)
      write(6,99) 'PDFmember',PDFmember
      write(6,*) '*******************************************'
      write(6,*)

c      write(6,*) '+ Name = ','PDFsets/'//PDFname
c      call InitPDFset('PDFsets/'//PDFname)
c      call InitPDFset(_PDFSETSPATH_//PDFname)
       call InitPDFsetByName(PDFname)
c      call InitPDFset(PDFname)
c      write(6,*) '+ PDF set succesfully initialized'
      
      if (PDFmember .lt. 0) then
        PDFerrors=.true.
        call numberPDF(maxPDFsets)
        if (maxPDFsets .gt. 1000) then
          write(6,*) 'ERROR: Max. number of error sets is 1000!'
          stop
        endif
        write(6,*)
        write(6,*) '****************************************'        
        write(6,*) '*        Calculating errors using      *'
        write(6,97) maxPDFsets
        write(6,*) '****************************************'
        call InitPDF(0)
        amz=alphasPDF(zmass)
        currentPDF=0
      else  
        call InitPDF(PDFmember)
        amz=alphasPDF(zmass)
      endif

c--- fill MCFM global variable "nlooprun" with the correct value
      call GetOrderAs(iorder)
      nlooprun=iorder+1
      
c--- rename pdlabel to get sensible output name
      pdlabel=PDFname(1:7)

      return
 
   97 format(' *        ',i4,' sets of error PDFs       *')
   98 format(' *   ',a7,' ',a29,' *')
   99 format(' *  ',a10,i3,'                          *')

  999 write(6,*) 'Error reading lhapdf.DAT'
      call flush(6)
      stop

      end
 

