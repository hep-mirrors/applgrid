      function GetNxx()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( MXX = 410 )
      PARAMETER ( MQ2 =  120 )
      COMMON/QCGRID/
     +     SCAX0,SCAQ0,XMICUT,QMICUT,QMACUT,RS2CUT,QMINAS,
     +     XXTAB(MXX),Q2TAB(MQ2),XHTAB(MXX),THRS34,THRS45,
     +     NXX,NQ2,NGRVER,IHTAB(MXX),NFMAP(MQ2),IQF2C(MQ2),
     +     IQF2B(MQ2),IQFLC(MQ2),IQFLB(MQ2),IFAILC(MXX,MQ2)
      integer getnxx
      getnxx = nxx
      return
      end

      function GetNq2()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( MXX = 410 )
      PARAMETER ( MQ2 =  120 )
      COMMON/QCGRID/
     +     SCAX0,SCAQ0,XMICUT,QMICUT,QMACUT,RS2CUT,QMINAS,
     +     XXTAB(MXX),Q2TAB(MQ2),XHTAB(MXX),THRS34,THRS45,
     +     NXX,NQ2,NGRVER,IHTAB(MXX),NFMAP(MQ2),IQF2C(MQ2),
     +     IQF2B(MQ2),IQFLC(MQ2),IQFLB(MQ2),IFAILC(MXX,MQ2)
      integer getnq2
      getnq2 = nq2
      return
      end

      function GetNgrver()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( MXX = 410 )
      PARAMETER ( MQ2 =  120 )
      COMMON/QCGRID/
     +     SCAX0,SCAQ0,XMICUT,QMICUT,QMACUT,RS2CUT,QMINAS,
     +     XXTAB(MXX),Q2TAB(MQ2),XHTAB(MXX),THRS34,THRS45,
     +     NXX,NQ2,NGRVER,IHTAB(MXX),NFMAP(MQ2),IQF2C(MQ2),
     +     IQF2B(MQ2),IQFLC(MQ2),IQFLB(MQ2),IFAILC(MXX,MQ2)
      integer GetNgrver
      GetNgrver = ngrver
      return
      end
      
      function GetXXtab(i)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( MXX = 410 )
      PARAMETER ( MQ2 =  120 )
      COMMON/QCGRID/
     +     SCAX0,SCAQ0,XMICUT,QMICUT,QMACUT,RS2CUT,QMINAS,
     +     XXTAB(MXX),Q2TAB(MQ2),XHTAB(MXX),THRS34,THRS45,
     +     NXX,NQ2,NGRVER,IHTAB(MXX),NFMAP(MQ2),IQF2C(MQ2),
     +     IQF2B(MQ2),IQFLC(MQ2),IQFLB(MQ2),IFAILC(MXX,MQ2)
      double precision Getxxtab
      Getxxtab = xxtab(i)
      return
      end
      
      function Getq2tab(i)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( MXX = 410 )
      PARAMETER ( MQ2 =  120 )
      COMMON/QCGRID/
     +     SCAX0,SCAQ0,XMICUT,QMICUT,QMACUT,RS2CUT,QMINAS,
     +     XXTAB(MXX),Q2TAB(MQ2),XHTAB(MXX),THRS34,THRS45,
     +     NXX,NQ2,NGRVER,IHTAB(MXX),NFMAP(MQ2),IQF2C(MQ2),
     +     IQF2B(MQ2),IQFLC(MQ2),IQFLB(MQ2),IFAILC(MXX,MQ2)
      double precision Getq2tab
      Getq2tab = q2tab(i)
      return
      end 


      function Getnfmap(i)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( MXX = 410 )
      PARAMETER ( MQ2 =  120 )
      COMMON/QCGRID/
     +     SCAX0,SCAQ0,XMICUT,QMICUT,QMACUT,RS2CUT,QMINAS,
     +     XXTAB(MXX),Q2TAB(MQ2),XHTAB(MXX),THRS34,THRS45,
     +     NXX,NQ2,NGRVER,IHTAB(MXX),NFMAP(MQ2),IQF2C(MQ2),
     +     IQF2B(MQ2),IQFLC(MQ2),IQFLB(MQ2),IFAILC(MXX,MQ2)
      integer Getnfmap
      Getnfmap = nfmap(i)
      return
      end



      function GetXMICUT()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( MXX = 410 )
      PARAMETER ( MQ2 =  120 )
      COMMON/QCGRID/
     +     SCAX0,SCAQ0,XMICUT,QMICUT,QMACUT,RS2CUT,QMINAS,
     +     XXTAB(MXX),Q2TAB(MQ2),XHTAB(MXX),THRS34,THRS45,
     +     NXX,NQ2,NGRVER,IHTAB(MXX),NFMAP(MQ2),IQF2C(MQ2),
     +     IQF2B(MQ2),IQFLC(MQ2),IQFLB(MQ2),IFAILC(MXX,MQ2)
      double precision GetXMICUT
      GetXMICUT = XMICUT
      return
      end

      function GetqMICUT()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( MXX = 410 )
      PARAMETER ( MQ2 =  120 )
      COMMON/QCGRID/
     +     SCAX0,SCAQ0,XMICUT,QMICUT,QMACUT,RS2CUT,QMINAS,
     +     XXTAB(MXX),Q2TAB(MQ2),XHTAB(MXX),THRS34,THRS45,
     +     NXX,NQ2,NGRVER,IHTAB(MXX),NFMAP(MQ2),IQF2C(MQ2),
     +     IQF2B(MQ2),IQFLC(MQ2),IQFLB(MQ2),IFAILC(MXX,MQ2)
      double precision GetqMICUT
      GetqMICUT = qMICUT
      return
      end
 
      function GetqMaCUT()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( MXX = 410 )
      PARAMETER ( MQ2 =  120 )
      COMMON/QCGRID/
     +     SCAX0,SCAQ0,XMICUT,QMICUT,QMACUT,RS2CUT,QMINAS,
     +     XXTAB(MXX),Q2TAB(MQ2),XHTAB(MXX),THRS34,THRS45,
     +     NXX,NQ2,NGRVER,IHTAB(MXX),NFMAP(MQ2),IQF2C(MQ2),
     +     IQF2B(MQ2),IQFLC(MQ2),IQFLB(MQ2),IFAILC(MXX,MQ2)
      double precision GetqMaCUT
      GetqMaCUT = qMaCUT
      return
      end

      function GetRS2CUT()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( MXX = 410 )
      PARAMETER ( MQ2 =  120 )
      COMMON/QCGRID/
     +     SCAX0,SCAQ0,XMICUT,QMICUT,QMACUT,RS2CUT,QMINAS,
     +     XXTAB(MXX),Q2TAB(MQ2),XHTAB(MXX),THRS34,THRS45,
     +     NXX,NQ2,NGRVER,IHTAB(MXX),NFMAP(MQ2),IQF2C(MQ2),
     +     IQF2B(MQ2),IQFLC(MQ2),IQFLB(MQ2),IFAILC(MXX,MQ2)
      double precision GetRS2CUT
      GetRS2CUT = RS2CUT
      return
      end

      function GetQMINAS()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( MXX = 410 )
      PARAMETER ( MQ2 = 120 )
      COMMON/QCGRID/
     +     SCAX0,SCAQ0,XMICUT,QMICUT,QMACUT,RS2CUT,QMINAS,
     +     XXTAB(MXX),Q2TAB(MQ2),XHTAB(MXX),THRS34,THRS45,
     +     NXX,NQ2,NGRVER,IHTAB(MXX),NFMAP(MQ2),IQF2C(MQ2),
     +     IQF2B(MQ2),IQFLC(MQ2),IQFLB(MQ2),IFAILC(MXX,MQ2)
      double precision GetQMINAS
      GetQMINAS = QMINAS
      return
      end

      


      function Getpdfqcd(i,ix,iq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NDFMAX = 20)
      PARAMETER ( MXX = 410 )
      PARAMETER ( MQ2 =  120 )
      COMMON/QCPASS/
     +     ALPHA0, Q0ALFA, ASLAST, QALAST,
     +     ALFASQ(MQ2), ALFAPQ(MQ2), ALFA2Q(MQ2),
     +     DELUP(MQ2), DELDN(MQ2), PDFQCD(MXX,MQ2,0:10),
     +     FNSQCD(MXX,MQ2),DNSQCD(MXX,MQ2),
     +     FSIQCD(MXX,MQ2),DSIQCD(MXX,MQ2),
     +     FGLQCD(MXX,MQ2),DGGQCD(MXX,MQ2),
     +     FSTORE(MXX,MQ2,31:30+NDFMAX),IDFAST(7,30),NDFAST,
     +     MARKFF(MXX,MQ2),MARKFH(MXX,MQ2),MARKQQ(MQ2),
     +     ISTFID(31:30+NDFMAX),IPDFID(31:30+NDFMAX),IEALFA(MQ2),
     +     IQL_LAST(10),IQ0_LAST(10),IQH_LAST(10)
      
      double precision Getpdfqcd
      Getpdfqcd = PDFQCD(ix,iq,i)
      return
      end 
      
      function Getpwgt(i,id,nf)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision Getpwgt
      COMMON /QCPWGT/ PWGT(0:10,0:30,3:5)
      Getpwgt = PWGT(i,id,nf)
      return
      end 
   

      subroutine opengridfile(gridname)
      implicit none
      character*(*) gridname
      integer qnerr

      qnerr=-1
      open(unit=2,status='old',file=gridname,
     .     form='unformatted',err=1)
      call QNREAD(2,1,qnerr)
 1    close(2)

      if (qnerr.ne.0) then
         write(*,*) 'Grid file problem: ',gridname
         if (qnerr.lt.0) then 
            write(*,*) 'Grid file does not exist'
            write(*,*) 'Calculating and creating grid file'
            call qnfilw(0,0)
            open(unit=2,status='unknown',file=gridname,
     .           form='unformatted')
            call QNDUMP(2)
            close(2)
         else
            write(*,*) 'Existing grid file is inconsistent'
            if (qnerr.eq.1)
     .           write(*,*) 'Defined grid different'
            if (qnerr.eq.2)
     .           write(*,*) 'Heavy quark weight table different'
            if (qnerr.eq.3)
     .           write(*,*) 'Charm mass different'
            if (qnerr.eq.4)
     .           write(*,*) 'Bottom mass different'
            stop
         endif
      endif
      end
      
      subroutine EVLSEA(name,IQ0,IQC,IQB,NQGRI)
      implicit none
      CHARACTER*(*) name
      integer iq0,iqc,iqb,nqgri
      real*8 f34,f45,f43,f54
      parameter(f34=1.d0/12.d0,f45=1.d0/20.d0,f43=-1.d0/12.d0,
     >          f54=-1.d0/20.d0)


      If(IQ0.le.IQC)then
         CAll EVPLUS(name,IQ0,1,IQC)
         CALL QADDSI(name,IQC,f34)
         CAll EVPLUS(name,IQC,IQC,IQB)
         CALL QADDSI(name,IQB,f45)
         CAll EVPLUS(name,IQB,IQB,NQGRI)
      else if(IQ0.gt.IQC.and.IQ0.le.IQB)then
         CAll EVPLUS(name,IQ0,IQC,IQB)
         CALL QADDSI(name,IQC,f43)
         CAll EVPLUS(name,IQC,1,IQC)
         CALL QADDSI(name,IQC,f34)
         CALL QADDSI(name,IQB,f45)
         CAll EVPLUS(name,IQB,IQB,NQGRI)
      else if(IQ0.gt.IQB)then
         CAll EVPLUS(name,IQ0,IQB,NQGRI)
         CALL QADDSI(name,IQB,f54)
         CAll EVPLUS(name,IQB,IQC,IQB)
         CALL QADDSI(name,IQB,f45)
         CALL QADDSI(name,IQC,f43)
         CAll EVPLUS(name,IQC,1,IQC)
         CALL QADDSI(name,IQC,f34)
      end if
      end


      subroutine EVOLCP(name,IQC,IQB,NQGRI)
      implicit none      
      CHARACTER*(*) name
      integer iqc,iqb,nqgri
      real*8 f4,f45
      parameter(f4=-1.d0/4.d0,f45=1.d0/20.d0)

c     First set to zero to avoid adding -1/4Singl at each iteration
      CAll QNPNUL(name)
      CALL QADDSI(name,IQC,f4)      
      CAll EVPLUS(name,IQC,IQC,IQB)
      CALL QADDSI(name,IQB,f45)
      CAll EVPLUS(name,IQB,IQB,NQGRI)
      end

      subroutine EVOLBP(name,IQB,NQGRI)
      implicit none
      CHARACTER*(*) name
      integer iqb,nqgri
      real*8 f5
      parameter(f5=-1.d0/5.d0)

c     First set to zero to avoid adding -1/5Singl at each iteration
      CAll QNPNUL(name)
      CALL QADDSI(name,IQB,f5)      
      CAll EVPLUS(name,IQB,IQB,NQGRI)
      end

