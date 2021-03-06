*********************************************************
*AUTHOR: FABIO MALTONI                                  *
*DATE  : 2/15/2004                                      *
*NOTES : PROGRAM GENERATED BY ToFortran-qq.m            *
*        AMPLITUDES CALCULATED BY ALBERTO FRIZZO        *
*********************************************************

      DOUBLE COMPLEX FUNCTION  A2q3g_mpppp(I1,I2,I3,I4,I5,za,zb)                        
* ---------------------------------------------------------------------
*                            1+ 2- 3+ 4+ 5+                            
* ---------------------------------------------------------------------


      IMPLICIT NONE
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      INTEGER I1,I2,I3,I4,I5
c      DOUBLE PRECISION       S(5,5)      
c      DOUBLE COMPLEX         ZA(5,5),ZB(5,5)
c      COMMON/PRODS/S,    ZA,     ZB     


       A2q3g_mpppp=(0d0,0d0)
       A2q3g_mpppp=A2q3g_mpppp+
     .        (za(i1,i2)**3*zb(i1,i3)**2*zb(i1,i4)*zb(i1,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i4)*za(i2,i5)) - 
     -  (za(i1,i2)**3*zb(i1,i3)*zb(i1,i4)**2*zb(i1,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i3)*za(i2,i5)) - 
     -  (za(i1,i2)**3*zb(i1,i3)*zb(i1,i4)**2*zb(i1,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i3)*za(i2,i5)) - 
     -  (za(i1,i2)**3*zb(i1,i3)*zb(i1,i4)**2*zb(i1,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i3)*za(i2,i5)) + 
     -  (za(i1,i2)**3*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)*za(i2,i4)) + 
     -  (za(i1,i2)**3*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)**2)/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)*za(i2,i4)) + 
     -  (za(i1,i2)**3*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)**2)/
     -   (s(i1,i3)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)*za(i2,i4)) + 
     -  (za(i1,i2)**2*za(i2,i3)*zb(i1,i3)**2*zb(i1,i5)*zb(i3,i4))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*s(i3,i4)*
     -     za(i2,i4)*za(i2,i5)) + 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i3,i4))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*s(i3,i4)*
     -     za(i2,i5)) + (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*
     -     zb(i1,i5)*zb(i3,i4))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i5)) + 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i3,i4))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i5)) + 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i3,i4))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i5)) + 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i3,i4))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i5)) + 
     -  (za(i1,i2)**2*za(i2,i4)*zb(i1,i4)**2*zb(i1,i5)*zb(i3,i4))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*s(i3,i4)*
     -     za(i2,i3)*za(i2,i5)) - 
     -  (za(i1,i2)**2*zb(i1,i3)*zb(i1,i5)**2*zb(i3,i4))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) - 
     -  (za(i1,i2)**2*zb(i1,i3)*zb(i1,i5)**2*zb(i3,i4))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) - 
     -  (za(i1,i2)**2*zb(i1,i3)*zb(i1,i5)**2*zb(i3,i4))/
     -   (s(i1,i2)*s(i3,i4)*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) - 
     -  (za(i1,i2)**2*zb(i1,i3)*zb(i1,i5)**2*zb(i3,i4))/
     -   (s(i1,i3)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) - 
     -  (za(i1,i2)**2*zb(i1,i3)*zb(i1,i5)**2*zb(i3,i4))/
     -   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) - 
     -  (za(i1,i2)**2*zb(i1,i4)*zb(i1,i5)**2*zb(i3,i4))/
     -   (s(i1,i2)*s(i3,i4)*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3))
       A2q3g_mpppp=A2q3g_mpppp
     .        -((za(i1,i2)**2*zb(i1,i4)*zb(i1,i5)**2*zb(i3,i4))/
     -     (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -       (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + 
     -         s(i2,i4) + s(i3,i4))*za(i2,i3))) - 
     -  (za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i5)*zb(i3,i4)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i5)) - 
     -  (za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i5)*zb(i3,i4)**2)/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i5)) - 
     -  (za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i5)*zb(i3,i4)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i5)) + 
     -  (za(i1,i2)*za(i2,i4)*zb(i1,i4)*zb(i1,i5)*zb(i3,i4)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i5)) + 
     -  (za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)**2*zb(i3,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i5)) + 
     -  (za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)**2*zb(i3,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i5)) - 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i3,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) - 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i3,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) - 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i3,i5))/
     -   (s(i1,i3)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) + 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i3,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i4)) - 
     -  (za(i1,i2)**2*zb(i1,i4)**2*zb(i1,i5)*zb(i3,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i3)) + 
     -  (za(i1,i2)*za(i2,i3)**2*zb(i1,i3)**2*zb(i3,i4)*zb(i3,i5))/
     -   (s(i1,i2)*s(i3,i4)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i4)*za(i2,i5)) - 
     -  (2*za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i4)*zb(i3,i4)*
     -     zb(i3,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i5)) - 
     -  (2*za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i4)*zb(i3,i4)*
     -     zb(i3,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i5)) + 
     -  (2*za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i4)*zb(i3,i4)*
     -     zb(i3,i5))/
     -   (s(i1,i2)*s(i3,i4)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i5)) + (za(i1,i2)*za(i2,i4)*zb(i1,i4)**2*zb(i3,i4)*
     -     zb(i3,i5))/
     -   (s(i1,i2)*s(i3,i4)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i5)) + (2*za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i3,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*s(i3,i4)*
     -     za(i2,i4)) + (2*za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i3,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) + 
     -  (2*za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i5)*zb(i3,i4)*
     -     zb(i3,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4))
       A2q3g_mpppp=A2q3g_mpppp+
     .        (2*za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i5)*zb(i3,i4)*
     -     zb(i3,i5))/
     -   (s(i1,i2)*s(i3,i4)*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) + 
     -  (2*za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i5)*zb(i3,i4)*
     -     zb(i3,i5))/
     -   (s(i1,i3)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) + 
     -  (2*za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i5)*zb(i3,i4)*
     -     zb(i3,i5))/
     -   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) + 
     -  (2*za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i5)*zb(i3,i4)*
     -     zb(i3,i5))/
     -   (s(i1,i2)*s(i3,i4)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i4)) + (2*za(i1,i2)*zb(i1,i4)*zb(i1,i5)*zb(i3,i4)*
     -     zb(i3,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*s(i3,i4)) + 
     -  (2*za(i1,i2)*zb(i1,i4)*zb(i1,i5)*zb(i3,i4)*zb(i3,i5))/
     -   (s(i1,i2)*s(i3,i4)*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) + (2*za(i1,i2)*zb(i1,i4)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i3,i5))/
     -   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) + (2*za(i1,i2)*zb(i1,i4)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i3,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))) + (2*za(i1,i2)*zb(i1,i4)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i3,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))) + (2*za(i1,i2)*zb(i1,i4)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i3,i5))/
     -   (s(i1,i2)*s(i3,i4)*(s(i3,i4) + s(i3,i5) + s(i4,i5))) + 
     -  (za(i1,i2)*za(i2,i5)*zb(i1,i5)**2*zb(i3,i4)*zb(i3,i5))/
     -   (s(i1,i2)*s(i3,i4)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i4)) + (za(i2,i3)**2*zb(i1,i3)*zb(i3,i4)**2*
     -     zb(i3,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i5)) + 
     -  (za(i2,i3)**2*zb(i1,i3)*zb(i3,i4)**2*zb(i3,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i5)) - 
     -  (za(i2,i3)*zb(i1,i5)*zb(i3,i4)**2*zb(i3,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))) + (za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i4)*
     -     zb(i3,i5)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) + 
     -  (za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i4)*zb(i3,i5)**2)/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) + 
     -  (za(i1,i2)*za(i2,i3)*zb(i1,i3)*zb(i1,i4)*zb(i3,i5)**2)/
     -   (s(i1,i3)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) + 
     -  (za(i1,i2)*za(i2,i5)*zb(i1,i4)*zb(i1,i5)*zb(i3,i5)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i4)) - 
     -  (za(i2,i3)**2*zb(i1,i3)*zb(i3,i4)*zb(i3,i5)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) - 
     -  (za(i2,i3)**2*zb(i1,i3)*zb(i3,i4)*zb(i3,i5)**2)/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4))
       A2q3g_mpppp=A2q3g_mpppp
     .        -((za(i2,i3)**2*zb(i1,i3)*zb(i3,i4)*zb(i3,i5)**2)/
     -     (s(i1,i2)*s(i3,i4)*
     -       (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + 
     -         s(i2,i4) + s(i3,i4))*za(i2,i4))) - 
     -  (za(i2,i3)**2*zb(i1,i3)*zb(i3,i4)*zb(i3,i5)**2)/
     -   (s(i1,i3)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) - 
     -  (za(i2,i3)**2*zb(i1,i3)*zb(i3,i4)*zb(i3,i5)**2)/
     -   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i4)) - 
     -  (za(i2,i3)*zb(i1,i4)*zb(i3,i4)*zb(i3,i5)**2)/
     -   (s(i1,i2)*s(i3,i4)*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) - (za(i2,i3)*zb(i1,i4)*zb(i3,i4)*
     -     zb(i3,i5)**2)/
     -   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) + (za(i2,i3)*za(i2,i5)*zb(i1,i5)*zb(i3,i4)*
     -     zb(i3,i5)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*s(i3,i4)*
     -     za(i2,i4)) + (za(i1,i2)**2*zb(i1,i3)**2*zb(i1,i4)*
     -     zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i5)) - 
     -  (za(i1,i2)**2*za(i2,i4)*zb(i1,i3)*zb(i1,i4)**2*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -     za(i2,i3)*za(i2,i5)) - 
     -  (za(i1,i2)**2*za(i2,i4)*zb(i1,i3)*zb(i1,i4)**2*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -     za(i2,i3)*za(i2,i5)) + 
     -  (za(i1,i2)**2*zb(i1,i3)**2*zb(i1,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i4)) + 
     -  (za(i1,i2)**2*zb(i1,i3)**2*zb(i1,i5)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i4)) - 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)) - 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)) - 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)) - 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i3)) - 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i3)) - 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i3)) - 
     -  (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*zb(i1,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -     za(i2,i3)) - (2*za(i1,i2)**2*zb(i1,i3)*zb(i1,i4)*
     -     zb(i1,i5)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -     za(i2,i3)) - (za(i1,i2)**2*za(i2,i5)*zb(i1,i3)*
     -     zb(i1,i5)**2*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -     za(i2,i3)*za(i2,i4))
       A2q3g_mpppp=A2q3g_mpppp
     .      -((za(i1,i2)**2*za(i2,i5)*zb(i1,i3)*zb(i1,i5)**2*zb(i4,i5))/
     -     (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -       za(i2,i3)*za(i2,i4))) + 
     -  (za(i1,i2)*za(i2,i3)*zb(i1,i3)**2*zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*s(i3,i4)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i5)) + (za(i1,i2)*za(i2,i3)*zb(i1,i3)**2*zb(i3,i4)*
     -     zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i5)) + (2*za(i1,i2)*za(i2,i4)*zb(i1,i3)*zb(i1,i4)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -     za(i2,i5)) + (2*za(i1,i2)*za(i2,i4)*zb(i1,i3)*zb(i1,i4)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -     za(i2,i5)) + (2*za(i1,i2)*za(i2,i4)*zb(i1,i3)*zb(i1,i4)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i5)) + 
     -  (2*za(i1,i2)*za(i2,i4)*zb(i1,i3)*zb(i1,i4)*zb(i3,i4)*
     -     zb(i4,i5))/
     -   (s(i1,i2)*s(i3,i4)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i5)) + (2*za(i1,i2)*za(i2,i4)*zb(i1,i3)*zb(i1,i4)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i5)) + (za(i1,i2)*za(i2,i4)**2*zb(i1,i4)**2*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*s(i3,i4)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i3)*za(i2,i5)) + 
     -  (za(i1,i2)*za(i2,i4)**2*zb(i1,i4)**2*zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i3)*za(i2,i5)) + 
     -  (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*s(i3,i4)) + 
     -  (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*s(i3,i4)*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5))
       A2q3g_mpppp=A2q3g_mpppp+
     .        (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)) + 
     -  (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*s(i3,i4)*(s(i3,i4) + s(i3,i5) + s(i4,i5))) + 
     -  (2*za(i1,i2)*zb(i1,i3)*zb(i1,i5)*zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*(s(i3,i4) + s(i3,i5) + s(i4,i5))) + 
     -  (2*za(i1,i2)*za(i2,i4)*zb(i1,i4)*zb(i1,i5)*zb(i3,i4)*
     -     zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*s(i3,i4)*
     -     za(i2,i3)) + (2*za(i1,i2)*za(i2,i4)*zb(i1,i4)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*s(i3,i4)*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)) + 
     -  (2*za(i1,i2)*za(i2,i4)*zb(i1,i4)*zb(i1,i5)*zb(i3,i4)*
     -     zb(i4,i5))/
     -   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)) + 
     -  (2*za(i1,i2)*za(i2,i4)*zb(i1,i4)*zb(i1,i5)*zb(i3,i4)*
     -     zb(i4,i5))/
     -   (s(i1,i2)*s(i3,i4)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i3)) + (2*za(i1,i2)*za(i2,i4)*zb(i1,i4)*zb(i1,i5)*
     -     zb(i3,i4)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i3)) + (za(i1,i2)*za(i2,i5)*zb(i1,i5)**2*zb(i3,i4)*
     -     zb(i4,i5))/
     -   (s(i1,i2)*s(i3,i4)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i3)) + (za(i1,i2)*za(i2,i5)*zb(i1,i5)**2*zb(i3,i4)*
     -     zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i3)) - (za(i2,i3)*za(i2,i4)*zb(i1,i3)*zb(i3,i4)**2*
     -     zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -     za(i2,i5)) - (za(i2,i3)*za(i2,i4)*zb(i1,i3)*zb(i3,i4)**2*
     -     zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -     za(i2,i5)) + (za(i2,i4)**2*zb(i1,i4)*zb(i3,i4)**2*
     -     zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i5)) + 
     -  (za(i2,i4)*zb(i1,i5)*zb(i3,i4)**2*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))) + (za(i2,i4)*zb(i1,i5)*zb(i3,i4)**2*
     -     zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))) + (za(i1,i2)*za(i2,i3)*zb(i1,i3)**2*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i4)) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i4)*zb(i3,i5)*
     -     zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4)))
       A2q3g_mpppp=A2q3g_mpppp+
     .        (2*za(i1,i2)*zb(i1,i3)*zb(i1,i4)*zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)) + 
     -  (2*za(i1,i2)*zb(i1,i3)*zb(i1,i4)*zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)) + 
     -  (2*za(i1,i2)*zb(i1,i3)*zb(i1,i4)*zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))) + (2*za(i1,i2)*zb(i1,i3)*zb(i1,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*(s(i3,i4) + s(i3,i5) + s(i4,i5))) + 
     -  (za(i1,i2)*za(i2,i4)*zb(i1,i4)**2*zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i3)) + (2*za(i1,i2)*za(i2,i5)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -     za(i2,i4)) + (2*za(i1,i2)*za(i2,i5)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -     za(i2,i4)) + (2*za(i1,i2)*za(i2,i5)*zb(i1,i3)*zb(i1,i5)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i4)) + 
     -  (2*za(i1,i2)*za(i2,i5)*zb(i1,i3)*zb(i1,i5)*zb(i3,i5)*
     -     zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i4)) + 
     -  (2*za(i1,i2)*za(i2,i5)*zb(i1,i3)*zb(i1,i5)*zb(i3,i5)*
     -     zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i4)) - (2*za(i1,i2)*za(i2,i5)*zb(i1,i4)*zb(i1,i5)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i3)) + 
     -  (2*za(i1,i2)*za(i2,i5)*zb(i1,i4)*zb(i1,i5)*zb(i3,i5)*
     -     zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i3)) + (za(i1,i2)*za(i2,i5)**2*zb(i1,i5)**2*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*(s(i3,i4) + s(i3,i5) + s(i4,i5))*
     -     za(i2,i3)*za(i2,i4)) - 
     -  (2*za(i2,i3)*zb(i1,i3)*zb(i3,i4)*zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) - (2*za(i2,i3)*zb(i1,i3)*zb(i3,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) - (2*za(i2,i3)*zb(i1,i3)*zb(i3,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*s(i3,i4)*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) - (2*za(i2,i3)*zb(i1,i3)*zb(i3,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4)))
       A2q3g_mpppp=A2q3g_mpppp+
     .        (-2*za(i2,i3)*zb(i1,i3)*zb(i3,i4)*zb(i3,i5)*zb(i4,i5))/
     -   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) - (2*za(i2,i3)*zb(i1,i3)*zb(i3,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))) - (2*za(i2,i3)*zb(i1,i3)*zb(i3,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))) - (2*za(i2,i3)*zb(i1,i3)*zb(i3,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)) - 
     -  (2*za(i2,i3)*zb(i1,i3)*zb(i3,i4)*zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)) - 
     -  (2*za(i2,i4)*zb(i1,i4)*zb(i3,i4)*zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*s(i3,i4)*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) - (2*za(i2,i4)*zb(i1,i4)*zb(i3,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) + (2*za(i2,i4)*zb(i1,i4)*zb(i3,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))) + (2*za(i2,i5)*zb(i1,i5)*zb(i3,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*s(i3,i4)) + 
     -  (2*za(i2,i5)*zb(i1,i5)*zb(i3,i4)*zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))) + (2*za(i2,i5)*zb(i1,i5)*zb(i3,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))) + (2*za(i2,i5)*zb(i1,i5)*zb(i3,i4)*
     -     zb(i3,i5)*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))) - (za(i2,i3)*za(i2,i5)*zb(i1,i3)*
     -     zb(i3,i5)**2*zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -     za(i2,i4)) - (za(i2,i3)*za(i2,i5)*zb(i1,i3)*zb(i3,i5)**2*
     -     zb(i4,i5))/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*s(i4,i5)*
     -     za(i2,i4)) + (za(i2,i5)*zb(i1,i4)*zb(i3,i5)**2*
     -     zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))) + (za(i2,i5)**2*zb(i1,i5)*zb(i3,i5)**2*
     -     zb(i4,i5))/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i4)) + 
     -  (za(i2,i5)**2*zb(i1,i5)*zb(i3,i5)**2*zb(i4,i5))/
     -   (s(i1,i2)*s(i4,i5)*
     -     (s(i1,i2) + s(i1,i4) + s(i1,i5) + s(i2,i4) + s(i2,i5) + 
     -       s(i4,i5))*za(i2,i4)) + 
     -  (za(i1,i2)*za(i2,i4)*zb(i1,i3)*zb(i1,i4)*zb(i4,i5)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)) + 
     -  (za(i1,i2)*za(i2,i4)*zb(i1,i3)*zb(i1,i4)*zb(i4,i5)**2)/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)) + 
     -  (za(i1,i2)*za(i2,i4)*zb(i1,i3)*zb(i1,i4)*zb(i4,i5)**2)/
     -   (s(i1,i3)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3))
       A2q3g_mpppp=A2q3g_mpppp
     .        -((za(i1,i2)*za(i2,i5)*zb(i1,i3)*zb(i1,i5)*zb(i4,i5)**2)/
     -     (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -       (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + 
     -         s(i2,i5) + s(i3,i5))*za(i2,i3))) - 
     -  (za(i1,i2)*za(i2,i5)*zb(i1,i3)*zb(i1,i5)*zb(i4,i5)**2)/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i3)) - 
     -  (za(i1,i2)*za(i2,i5)*zb(i1,i3)*zb(i1,i5)*zb(i4,i5)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i3)) - 
     -  (za(i2,i4)*zb(i1,i3)*zb(i3,i4)*zb(i4,i5)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) - (za(i2,i4)*zb(i1,i3)*zb(i3,i4)*
     -     zb(i4,i5)**2)/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) - (za(i2,i4)*zb(i1,i3)*zb(i3,i4)*
     -     zb(i4,i5)**2)/
     -   (s(i1,i2)*s(i3,i4)*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) - (za(i2,i4)*zb(i1,i3)*zb(i3,i4)*
     -     zb(i4,i5)**2)/
     -   (s(i1,i3)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) - (za(i2,i4)*zb(i1,i3)*zb(i3,i4)*
     -     zb(i4,i5)**2)/
     -   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))) - (za(i2,i4)**2*zb(i1,i4)*zb(i3,i4)*
     -     zb(i4,i5)**2)/
     -   (s(i1,i2)*s(i3,i4)*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)) - 
     -  (za(i2,i4)**2*zb(i1,i4)*zb(i3,i4)*zb(i4,i5)**2)/
     -   (s(i3,i4)*(s(i1,i3) + s(i1,i4) + s(i3,i4))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i4) + s(i2,i3) + s(i2,i4) + 
     -       s(i3,i4))*za(i2,i3)) + 
     -  (za(i2,i4)*za(i2,i5)*zb(i1,i5)*zb(i3,i4)*zb(i4,i5)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*s(i3,i4)*
     -     za(i2,i3)) + (za(i2,i5)*zb(i1,i3)*zb(i3,i5)*
     -     zb(i4,i5)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))) + (za(i2,i5)*zb(i1,i3)*zb(i3,i5)*
     -     zb(i4,i5)**2)/
     -   (s(i1,i3)*(s(i1,i2) + s(i1,i3) + s(i2,i3))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))) - (za(i2,i5)**2*zb(i1,i5)*zb(i3,i5)*
     -     zb(i4,i5)**2)/
     -   (s(i1,i2)*(s(i1,i2) + s(i1,i5) + s(i2,i5))*
     -     (s(i1,i2) + s(i1,i3) + s(i1,i5) + s(i2,i3) + s(i2,i5) + 
     -       s(i3,i5))*za(i2,i3))

      END


