      block data electroweak_input
************************************************************************
*     Calculational scheme for EW couplings                            *
************************************************************************
c
c     ewscheme=-1  : Old MCFM default 
c                    input values = Gf,alpha(m_Z),m_W,m_Z
c                    output values = sin^2(theta_W),mtop
c
c     ewscheme=0   : Old MadEvent default (= AlpGen with iewopt=2)
c                    input values = sin^2(theta_W),alpha(m_Z),m_Z
c                    output values = m_W,Gf.
c
c     ewscheme=1   : New MCFM default, also Madevent default, "G_mu scheme"
c                    = LUSIFER and AlpGen (iewopt=3) defaults
c                    input values = G_F,m_Z,m_W
c                    output values = sin^2(theta_W),alpha(m_Z).
c
c     ewscheme=2   : input  values = G_F,sin^2(theta_W),alpha(m_Z)
c                    output values = m_W,m_Z.
c
c     ewscheme=3   : User choice. All parameters are left as they are
c                    input here. You have to know what you're doing.
c
      implicit none
      include 'ewinput.f'
      data ewscheme  / +1                  /   ! Chooses EW scheme
      data Gf_inp    / 1.16639d-5          /   ! G_F
      data aemmz_inp / 7.7585538055706d-03 /   ! alpha_EM(m_Z)=1/128.89
      data xw_inp    / 0.2312d0            /   ! sin^2(theta_W)
      data wmass_inp / 80.419d0            /   ! W mass
      data zmass_inp / 91.188d0            /   ! Z mass
      end
************************************************************************


************************************************************************
*     Masses, widths and initial-state flavour information             *
************************************************************************
      block data block_properties
      implicit none
      include 'masses.f'
      include 'nflav.f'
      include 'nores.f'
c--- if true, nores removes all of the gg contribution
      data nores/.false./
c--- Masses: note that "mtausq", "mcsq" and "mbsq" are typically used
c--- throughout the program to calculate couplings that depend on the
c--- mass, while "mtau","mc" and "mb" are the masses that appear in
c--- the rest of the matrix elements and phase space (and may be set
c--- to zero in the program, depending on the process number) 
      data mtausq,mcsq,mbsq/3.157729d0,2.25d0,21.3444d0/
      data mtau/1.777d0/
      data mc,mb,mt/1.5d0,4.62d0,170.9d0/
c--- Widths: note that the top width is calculated in the program
      data wwidth,zwidth/2.06d0,2.49d0/
      data tauwidth/2.269d-12/
c--- Number of active flavours in the initial state: this parameter
c--- may be changed in the program for some processes
      data nflav/5/
c--- Masses below here are currently unused      
      data md,mu,ms/5d-3,5d-3,1d-1/
      data mel,mmu/0.510997d-3,0.105658389d0/
      end
************************************************************************


************************************************************************
*     CKM matrix entries                                               *
************************************************************************
      block data block_ckm
      implicit none
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb
      data  Vud  ,  Vus  ,  Vub  ,
     .      Vcd  ,  Vcs  ,  Vcb
     .   /0.975d0,0.222d0,0.000d0,
     .    0.222d0,0.975d0,0.000d0/
      end
************************************************************************


************************************************************************
*     Relevant for the H+b process only :                              *
*       mb_msbar: the value of the running b-mass, evaluated at the    * 
*                 pole mass. For negative values, calculated from mb   *
*       susycoup: the deviation of the Higgs coupling from the         *
*                 Standard Model value (S.M. = 1d0)                    *
************************************************************************
      block data block_bH
      implicit none
      include 'mb_msbar.f'
      include 'susycoup.f'
      data mb_msbar/4.25d0/
      data susycoup/1d0/
      end
************************************************************************


************************************************************************
*     Dim. Reg. parameter epsilon, used for checking the proper        *
*      operation of the NLO code in the program                        *
************************************************************************
      block data block_epinv
      implicit none
      include 'epinv.f'
      include 'epinv2.f'
      data epinv/ 1d3/
      data epinv2/1d3/
      end
************************************************************************

