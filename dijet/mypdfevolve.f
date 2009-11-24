       subroutine mypdfevolve(x,Q,f)
       implicit none
       real*8 x,Q,f(-6:6)
       integer i

C       logical debug
C       data debug/.false./ 

C       if (debug) write(6,*) '  x= ',x,' Q= ',Q

       call evolvePDF(x,Q,f)

C       if (debug) then
C        do i = -6,6
C          write(6,*) 'in: i= ',i,' f= ',f(i)
C        enddo
C       endif      

C        do i = -6,6
C         if (f(i).gt.100000) then
C          write(6,*) '**mypdf: x= ',x,' Q= ',Q,' f= ',f(i),' reset=0 '
C         endif
C        enddo

       if (x.eq.1.) then
        do i = -6,6
         f(i)=0.
        enddo
       endif

C       if (debug) then
C        do i = -6,6
C          write(6,*) 'out: i= ',i,' f= ',f(i)
C        enddo
C       endif      

       return
       end
