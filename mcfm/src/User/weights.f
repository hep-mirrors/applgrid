cc
cc    initializes contribution from different
cc    subprocesses to zero
cc
      subroutine iniContribution(contrib)
      implicit none
      include 'APPLinclude.f'
      double precision contrib(0:nSubProcess-1)
      integer iProcess

      do iProcess = 0, nSubProcess-1
         contrib(iProcess) = 0d0
      enddo

      return
      end
cc*****************************************************
cc
cc    scales contribution from different
cc    subprocesses by factor (phase space, branching, etc)
cc
      subroutine scaleContribution(factor,contrib)
      implicit none
      include 'APPLinclude.f'
      double precision factor
      double precision contrib(0:nSubProcess-1)
      integer iProcess

      do iProcess = 0, nSubProcess-1
         contrib(iProcess) = factor * contrib(iProcess)
      enddo

      return
      end
cc*****************************************************
cc
cc    fills weight to array contrib accordong 
cc    to the flavour combination 
cc
      subroutine collectWeights(iflav1,iflav2,currentWeight,contrib)
      implicit none
      include 'APPLinclude.f'
      integer iflav1, iflav2
      double precision contrib(0:nSubProcess-1)
      double precision currentWeight, factor
      integer iSubProcess
      logical debug222
c  default values
      debug222 = .false.
      factor = 1d0

      call makeDecision(iflav1,iflav2,iSubProcess,factor)

      if (iSubProcess.ge.0) then
         contrib(iSubProcess) = contrib(iSubProcess) + 
     &                          currentWeight*factor
         
         if(debug222.and.(currentWeight.ne.0d0))then
            print*,"*** collectWeights:   cw = ",currentWeight
            print*,"*** collectWeights:    f = ",factor
            print*,"*** collectWeights:    w = ",currentWeight*factor
            print*,"*** collectWeights: subP =  ",iSubProcess
         endif
         
      endif
      
      return
      end
cc*****************************************************
cc  defines the subprocess number for flavour combination (iflav1, iflav2)
cc  and calculates factor from CKM matrix element
cc
      subroutine makeDecision(iflav1,iflav2,iProcess,factor)
      implicit none
      include 'constants.f'
      include 'ckm.f'
      integer iflav1, iflav2, iProcess
      double precision factor
      integer nproc
      common /nproc/ nproc
c     default value (do not fill array)
      iProcess = -1
      factor = 1d0
c--->>>>>b bbar production process      
      if (nproc.eq.158)then
         if ((iflav1.eq.0).or.(iflav2.eq.0)) then
            iProcess = 0
         endif
         if ((iflav1.eq.1).and.(iflav2.eq.-1)) then
            iProcess = 1
         endif
      endif
c--->>>>> W_only production W+
c
c   0: dbar-quark <--> u-quark
c   1:    u-quark <--> dbar-quark
c   2: dbar-quark <--> Gluon
c   3:   u-quark  <--> Gluon
c   4:     Gluon  <--> dbar-quark
c   5:     Gluon  <--> u-quark
c
      if (nproc.eq.1)then
         
         if((iflav1.eq.0).and.(iflav2.eq. 2))then
            factor = 1d0/Vsum(iflav2)
            iProcess = 5
         endif
         if((iflav1.eq.0).and.(iflav2.eq.-1))then
            factor = 1d0/Vsum(iflav2)
            iProcess = 4
         endif
         if((iflav2.eq.0).and.(iflav1.eq. 2))then
            factor = 1d0/Vsum(iflav1)
            iProcess = 3
         endif
         if((iflav2.eq.0).and.(iflav1.eq.-1))then
            factor = 1d0/Vsum(iflav1)
            iProcess = 2
         endif
         if ((iflav1.eq.2).and.(iflav2.eq.-1)) then
            factor = 1d0/Vsq(iflav1,iflav2)
            iProcess = 1
         endif
         if ((iflav1.eq.-1).and.(iflav2.eq.2)) then
            factor = 1d0/Vsq(iflav1,iflav2)
            iProcess = 0
         endif
      endif

c---  >>>>> W_only production W-
c
c   0: d-quark    <-->  ubar-quark
c   1: ubar-quark <-->  d-quark
c   2: d-quark    <-->  Gluon
c   3: ubar-quark <-->  Gluon
c   4: Gluon      <-->  d-quark
c   5: Gluon      <-->  ubar-quark
c

      if (nproc.eq.6)then
         
         if((iflav1.eq.0).and.(iflav2.eq. -2))then
            factor = 1d0/Vsum(iflav2)
            iProcess = 5
         endif
         if((iflav1.eq.0).and.(iflav2.eq.1))then
            factor = 1d0/Vsum(iflav2)
            iProcess = 4
         endif
         if((iflav2.eq.0).and.(iflav1.eq.-2))then
            factor = 1d0/Vsum(iflav1)
            iProcess = 3
         endif
         if((iflav2.eq.0).and.(iflav1.eq.1))then
            factor = 1d0/Vsum(iflav1)
            iProcess = 2
         endif
         if ((iflav1.eq.-2).and.(iflav2.eq.1)) then
            factor = 1d0/Vsq(iflav1,iflav2)
            iProcess = 1
         endif
         if ((iflav1.eq.1).and.(iflav2.eq.-2)) then
            factor = 1d0/Vsq(iflav1,iflav2)
            iProcess = 0
         endif

      endif
c
c--->>>>> Z_only production
c
c   0:       UP <--> UPbar
c   1:     DOWN <--> DOWNbar
c   2:    UPbar <--> UP
c   3   DOWNbar <--> DOWN
c   4:    Gluon <--> UP
c   5:    Gluon <--> UPbar
c   6:    Gluon <--> DOWN
c   7:    Gluon <--> DOWNbar
c   8:       UP <--> Gluon
c   9:    UPbar <--> Gluon
c  10:     DOWN <--> Gluon
c  11:  DOWNbar <--> Gluon

      if ((nproc.eq.31)) then

         if (iflav2.eq.0) then         
cc
            if     (iflav1.eq.-1) then
               iProcess = 11
            elseif (iflav1.eq. 1) then
               iProcess = 10
            elseif (iflav1.eq.-2) then
               iProcess = 9
            elseif (iflav1.eq. 2) then
               iProcess = 8
            endif
ccc
         elseif (iflav1.eq.0) then         
cc
            if     (iflav2.eq.-1) then
               iProcess = 7
            elseif (iflav2.eq. 1) then
               iProcess = 6
            elseif (iflav2.eq.-2) then
               iProcess = 5
            elseif (iflav2.eq. 2) then
               iProcess = 4
            endif
ccc
         elseif ((iflav1*iflav2).ne.0) then

            if     (iflav1.eq.-1) then
               iProcess = 3
            elseif (iflav1.eq.-2) then
               iProcess = 2
            elseif (iflav1.eq. 1) then
               iProcess = 1
            elseif (iflav1.eq. 2) then
               iProcess = 0
            endif

         endif

         factor =1d0
      endif                     ! if( nproc == 31)
c     
      return
      end
c*****************************************************
c*****************************************************
c*****************************************************
c*********      REAL CONTRIBUTION    *****************
c*****************************************************
c*****************************************************
c*****************************************************
c    initializes contribution from different
c    subprocesses to zero for REAL CONTRIBUTION
c
      subroutine iniRContribution(index,contrib)
      implicit none
      include 'constants.f'
      include 'ptilde.f'
      include 'APPLinclude.f'
      integer index
      double precision contrib(0:maxd,0:nSubProcess-1)
      double precision currentContrib(0:nSubProcess-1)
      integer iProcess

      call iniContribution(currentContrib)
      
      do iProcess = 0, nSubProcess-1
         contrib(index,iProcess) = currentContrib(iProcess)
      enddo
      
      return
      end
cc*****************************************************
cc
cc    scales contribution from different
cc    subprocesses by factor (phase space, branching, etc)
cc    for REAL CONTRIBUTION
cc
      subroutine scaleRContribution(index,factor,contrib)
      implicit none
      include 'constants.f'
      include 'ptilde.f'
      include 'APPLinclude.f'
      integer index
      double precision factor
      double precision contrib(0:maxd,0:nSubProcess-1)
      double precision currentContrib(0:nSubProcess-1)
      integer iProcess

      do iProcess=0, nSubProcess-1
         currentContrib(iProcess) = contrib(index,iProcess)
      enddo
      
      call scaleContribution(factor,currentContrib)

      do iProcess=0, nSubProcess-1
         contrib(index,iProcess) = currentContrib(iProcess)
      enddo
      
      return
      end
cc*****************************************************
cc
cc    fills weight to array contrib accordong 
cc    to the flavour combination 
cc    for REAL CONTRIBUTION
cc
      subroutine collectRWeights(
     &     index,iflav1,iflav2,
     &     currentWeight,contrib)
      implicit none
      include 'constants.f'
      include 'ptilde.f'
      include 'APPLinclude.f'
      integer index
      integer iflav1, iflav2
      double precision currentWeight
      double precision contrib(0:maxd,0:nSubProcess-1)
      double precision tmpContrib(0:nSubProcess-1)
      double precision currentContrib(0:nSubProcess-1)
      integer iProcess
      
cc      print*,"get w = ",currentWeight," in ( ",iflav1," , ",iflav2," )"


      do iProcess=0, nSubProcess-1
         tmpContrib(iProcess) = contrib(index,iProcess)
      enddo

      call iniContribution(currentContrib)
      
      call collectWeights(iflav1,iflav2,currentWeight,currentContrib)
      
      do iProcess=0, nSubProcess-1
         contrib(index,iProcess) = 
     &        tmpContrib(iProcess) + currentContrib(iProcess)
      enddo
      
      return
      end
cc*****************************************************
cc*****************************************************
cc*****************************************************
cc*****************************************************
cc         NEEDED FOR SUBPROCESS REFERENCE HISTOGRAMS
cc*****************************************************
cc*****************************************************
cc*****************************************************
cc*****************************************************
cc*****************************************************
      integer function getProcess(iflav1,iflav2)
      implicit none
      integer iflav1, iflav2
      integer iProcess
      integer nproc
      common /nproc/ nproc
ccc   default value
      getProcess = -1
c------->>>> b bbar production process      

c------->>>> w production
c W+
      if (nproc.eq.1)then
         if((iflav1.eq.0).and.(iflav2.gt.0)) getProcess = 5
         if((iflav1.eq.0).and.(iflav2.lt.0)) getProcess = 4
         if((iflav2.eq.0).and.(iflav1.gt.0)) getProcess = 3
         if((iflav2.eq.0).and.(iflav1.lt.0)) getProcess = 2
         if((iflav1.gt.0).and.(iflav2.lt.0)) getProcess = 1
         if((iflav1.lt.0).and.(iflav2.gt.0)) getProcess = 0
      endif
c W-
      if (nproc.eq.6)then
         if((iflav1.eq.0).and.(iflav2.lt.0)) getProcess = 5
         if((iflav1.eq.0).and.(iflav2.gt.0)) getProcess = 4
         if((iflav2.eq.0).and.(iflav1.lt.0)) getProcess = 3
         if((iflav2.eq.0).and.(iflav1.gt.0)) getProcess = 2
         if((iflav1.lt.0).and.(iflav2.gt.0)) getProcess = 1
         if((iflav1.gt.0).and.(iflav2.lt.0)) getProcess = 0
      endif
c------->>>> z production
      if (nproc.eq.31) then
         if     (iflav2.eq.0) then
ccccccccc
            if (mod(abs(iflav1),2).ne.0) then
               if (iflav1.lt.0) then 
                  getProcess = 11
               else
                  getProcess = 10
               endif
            else
               if (iflav1.lt.0) then 
                  getProcess = 9
               else
                  getProcess = 8
               endif
            endif
ccccccccc
         elseif (iflav1.eq.0) then
            
            if (mod(abs(iflav2),2).ne.0) then
               if (iflav2.lt.0) then 
                  getProcess = 7
               else
                  getProcess = 6
               endif
            else
               if (iflav2.lt.0) then 
                  getProcess = 5
               else
                  getProcess = 4
               endif
            endif

         elseif ((iflav1*iflav2).ne.0) then
            if (mod(abs(iflav1),2).ne.0) then 
               if (iflav1.gt.0) then
                  getProcess = 1
               else
                  getProcess = 3
               endif
            else
               if (iflav1.gt.0) then
                  getProcess = 0
               else
                  getProcess = 2
               endif
            endif
         endif
      endif
      
      return
      end
cc
cc
      subroutine collectReference(iflav1,iflav2,curWeight,contrib)
      implicit none
      include 'APPLinclude.f'
      integer iflav1, iflav2
      double precision curWeight
      double precision contrib(0:nSubProcess-1)
      integer iProcess, getProcess
      logical debug222
c  default values
      debug222 = .false.
cc
      iProcess = getProcess(iflav1, iflav2)
      
      contrib(iProcess) = contrib(iProcess) + curWeight
      
      if(debug222.and.(curWeight.ne.0d0))then 
         print*,"      ProcRef-->> w: ",iflav1,iflav2," w = ",curWeight
         print*,"      ProcRef-->>  filled  ",curWeight," in ",iProcess
      endif
      return
      end
cc
cc
cc      
      subroutine collectRReference(index,iflav1,iflav2,
     &                             curWeight,contrib)
      implicit none
      include 'APPLinclude.f'
      include 'constants.f'
      include 'ptilde.f'
      integer index,iflav1, iflav2
      double precision curWeight
      double precision contrib(0:maxd,0:nSubProcess-1)
      integer iProcess, getProcess
      logical debug222
ccc  default values
      debug222 = .false.

      iProcess = getProcess(iflav1,iflav2)

      if(debug222)print*,"weights: ",iflav1,iflav2," w = ",curWeight
      
      contrib(index,iProcess) = contrib(index,iProcess) + curWeight
      
      if(debug222)then
         print*,"       filled",curWeight," in ",iProcess,
     &                " for dipole =", index
      endif
      
      return
      end

