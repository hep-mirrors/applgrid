module f77_pdftab
  use types; use consts_dp
  use pdf_tabulate_new
  use convolution; use pdf_general; use conv_objects
  use holders; use pdf_general; use dglap_choices
  implicit none

  !! holds information about the grid
  type(grid_def),     save :: gd, gdarray(3)

  !! holds the splitting functions
  type(sigma_holder), save :: sh

  !! 0 is main pdf table, while i=1:3 contain convolutions with the
  !! i-loop splitting function
  type(pdftab), save :: tables(0:3)
  logical,      save :: setup_done(0:3) = .false.
  integer,      save :: setup_nf(3)     = 0
end module f77_pdftab


!======================================================================
!! initialise the underlying grid, splitting functions and pdf-table
!! objects, using the dy and nloop parameters as explained below.
subroutine dglapStart(dy,nloop)
  use f77_pdftab
  implicit none
  real(dp), intent(in) :: dy     !! internal grid spacing: 0.1 is a sensible value
  integer,  intent(in) :: nloop  !! the maximum number of loops we'll want (<=3)
  !-------------------------------------
  real(dp)          :: ymax
  integer           :: order

  ! initialise our grids
  ! specify the maximum value of log(1/x)
  ymax = 12.0_dp
  ! the internal interpolation order (with a minus sign allows
  ! interpolation to take fake zero points beyond x=1 -- convolution
  ! times are unchanged, initialisation time is much reduced and
  ! accuracy is slightly reduced)
  order = -5 
  ! Now create a nested grid
  call conv_InitGridDef(gdarray(3),dy/9.0_dp,0.5_dp, order=order)
  call conv_InitGridDef(gdarray(2),dy/3.0_dp,2.0_dp, order=order)
  call conv_InitGridDef(gdarray(1),dy,       ymax  ,order=order)
  call conv_InitGridDef(gd,gdarray(1:3),locked=.true.)

  ! create the tables that will contain our copy of the user's pdf
  ! as well as the convolutions with the pdf.
  ! call pdftab_AllocTab(gd, tables(:), Qmin=1.0_dp, Qmax=10000.0_dp, & 
    call pdftab_AllocTab(gd, tables(:), Qmin=1.0_dp, Qmax=28000.0_dp, & 
 !& dlnlnQ = min(dy,0.1_dp), freeze_at_Qmin=.true.)
  & dlnlnQ = 0.03_dp, freeze_at_Qmin=.true.)

  ! initialise splitting-function holder
  call holder_InitSigma(gd,sh,factscheme=factscheme_MSbar,&
       &                      nloop=nloop,nflo=3,nfhi=6)
  ! choose a sensible default number of flavours.
  call holder_SetNf(sh,nflcl=5)

  ! indicate the pdfs and convolutions have not been initialised...
  setup_done = .false.
end subroutine dglapStart


!======================================================================
!! Given a pdf_subroutine with the interface shown below, initialise
!! our internal pdf table.
subroutine dglapAssign(pdf_subroutine)
  use f77_pdftab
  implicit none
  interface
     subroutine pdf_subroutine(x,Q,res)
       use types; implicit none
       real(dp), intent(in)  :: x,Q
       real(dp), intent(out) :: res(*)
     end subroutine pdf_subroutine
  end interface
  !-----------------------------------

  call pdftab_InitTab_LHAPDF(tables(0), pdf_subroutine)
  setup_done(0)  = .true.
  setup_done(1:) = .false.
end subroutine dglapAssign


!======================================================================
!! Return in f(-6:6) the value of the internally stored pdf at the
!! given x,Q, with the usual LHApdf meanings for the indices -6:6.
subroutine dglapEval(x,Q,f)
  use f77_pdftab
  implicit none
  real(dp), intent(in)  :: x, Q
  real(dp), intent(out) :: f(-6:6)
  
  call pdftab_ValTab_xQ(tables(0),x,Q,f)
end subroutine dglapEval



!======================================================================
!! Return in f(-6:6) the value of 
!!
!!    [P(iloop,nf) \otimes pdf] (x,Q)
!!
!! where P(iloop,nf) is the iloop-splitting function for the given
!! value of nf, and pdf is our internally stored pdf.
!!
!! The normalisation is such that the nloop dglap evolution equation is
!!
!!     dpdf/dlnQ^2 = sum_{iloop=1}^nloop 
!!                        (alphas/(2*pi))^iloop * P(iloop,nf) \otimes pdf
!!
!! Note that each time nf changes relative to a previous call for the
!! same iloop, the convolution has to be repeated for the whole
!! table. So for efficient results when requiring multiple nf values,
!! calls with the same nf value should be grouped together.
!!
!! In particular, for repeated calls with the same value of nf, the
!! convolutions are carried out only on the first call (i.e. once for
!! each value of iloop). Multiple calls with different values for
!! iloop can be carried out without problems.
!!
subroutine dglapEvalSplit(x,Q,iloop,nf,f)
  use f77_pdftab; use warnings_and_errors
  implicit none
  real(dp), intent(in)  :: x, Q
  integer,  intent(in)  :: iloop, nf
  real(dp), intent(out) :: f(-6:6)
  
  if (.not. setup_done(iloop) .or. setup_nf(iloop) /= nf) then
     if (iloop > size(sh%allP,dim=1) .or. iloop < 1) &
          &call wae_error('dglapeval_split','illegal value for iloop:',&
          &intval=iloop)

     if (nf < lbound(sh%allP,dim=2) .or. nf > ubound(sh%allP,dim=2)) &
          &call wae_error('dglapeval_split','illegal value for nf:',&
          &intval=nf)

     tables(iloop)%tab = sh%allP(iloop, nf) .conv. tables(0)%tab
     
     setup_done(iloop) = .true.
     setup_nf(iloop)   = nf
  end if
  
  call pdftab_ValTab_xQ(tables(iloop),x,Q,f)
end subroutine dglapEvalSplit


