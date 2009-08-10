!-----------------------------------------------------------
!! Module that will hopefully make tabulation easy
!!
!! Currently does not know anything about nf thresholds.
!!
!! But some form of support for these might be useful... Question is:
!! how can it be arranged in a fairly clean manner, without breaking
!! things that already use the system (e.g. caesar resum)?
!!
!! One possibility is to use an as_handle as a source of information.
!! One can then either take a copy of the as_handle and use its routines
!! or else just extract the relevant useful information
!!
!! Aim is to be able to write -- i.e. convolutions with full nf info
!!
!!   Ptab%tab = sh%P(iloop,tab%nf_int(:)) .conv. tab%tab
!!
!! This probably requires some modifications also to the convolution 
!! routines (which do not currently have this overloading for .conv.).
!!
!! [NB at some point maybe include nf association in some simpler manner?]
module pdf_tabulate_new
  use types; use consts_dp
  use convolution; use conv_objects
  use pdf_representation; use pdf_general
  use interpolation
  use warnings_and_errors
  !-- needed only for pdftab_InitTabEvolve [separate it out one day?]
  use as_server; use evolution; use holders
  implicit none
  private

  
  type pdfseginfo
     real(dp) :: lnlnQ_lo, lnlnQ_hi, dlnlnQ
     integer :: ilnlnQ_lo, ilnlnQ_hi
  end type pdfseginfo
  public :: pdfseginfo

  type pdftab
     ! basic elements of a pdftab, common regardless of whether we
     ! additionally have the nf segments...
     type(grid_def) :: gd
     real(dp) :: default_dlnlnQ
     real(dp) :: lnlnQ_min, lnlnQ_max, lambda_eff
     real(dp), pointer :: tab(:,:,:)
     real(dp), pointer :: lnlnQ_vals(:)
     integer :: nQ
     logical :: freeze_at_Qmin
     ! this is useful only in absence of nf info.
     real(dp) :: dlnlnQ
     !
     ! Stuff to do with variable nf and alpha_s; not always available.
     ! In cases with variable nf, the table will be broken into multiple
     ! segments, each one of which potentially different spacings.
     !
     logical :: nf_info_associated
     integer :: nflo, nfhi
     type(pdfseginfo), pointer :: seginfo(:)
     integer, pointer :: nf_int(:)
     real(dp), pointer :: as2pi(:)
     !
     ! Elements needed in case we want to do precalculation of
     ! of evolution. Not always available.
     type(evln_operator), pointer :: evops(:)
     integer                      :: StartScale_iQlo
     real(dp)                     :: StartScale
  end type pdftab
  public :: pdftab

  !-- for calculating ln ln Q/lambda_eff
  real(dp), parameter :: default_lambda_eff = 0.1_dp
  real(dp), parameter :: default_dlnlnQ = 0.1_dp
  real(dp), parameter :: warn_tolerance = 1e-3_dp
  ! used in various contexts for deciding when an interval is
  ! sufficiently small that it can be ignored...
  real(dp), parameter :: min_dlnlnQ_singleQ = 1e-10_dp
  !integer, parameter :: lnlnQ_order = 3
  !integer, parameter :: lnlnQ_order = 2
   integer, parameter :: lnlnQ_order = 4

  interface pdftab_AllocTab
     module procedure pdftab_AllocTab_, pdftab_AllocTab_fromorig,&
          & pdftab_AllocTab_1d, pdftab_AllocTab_fromorig_1d
  end interface

  interface pdftab_AssocNfInfo
     module procedure pdftab_AssocNfInfo, pdftab_AssocNfInfo_1d
  end interface

  interface pdftab_InitTabSub
     module procedure pdftab_InitTabSub_, pdftab_InitTabSub_iset
  end interface

  public :: pdftab_InitTab_LHAPDF

  interface pdftab_InitTabEvolve
     module procedure pdftab_InitTabEvolve_frompre, pdftab_InitTabEvolve
  end interface

  interface pdftab_DelTab
     module procedure pdftab_DelTab_0d, pdftab_DelTab_1d
  end interface

  public :: pdftab_AllocTab, pdftab_InitTabSub
  public :: pdftab_AssocNfInfo
  public :: pdftab_InitTabEvolve, pdftab_PreEvolve
  public :: pdftab_TabEvolveGen
  public :: pdftab_ValTab_yQ, pdftab_ValTab_xQ
  public :: pdftab_DelTab

contains

  !---------------------------------------------------------
  !! Allocate a pdftab, which covers a range Qmin to Qmax, using
  !! tabulation uniform in (ln ln Q/Lambda), where Lambda is currently
  !! a module parameter.
  !!
  !! If freeze_at_Qmin is present and .true., distributions are frozen
  !! at their Qmin value for Q < Qmin. Otherwise they are set to zero
  !! there.
  !!
  !! More sensible extrapolation beyond Q range offers scope for future 
  !! improvement here!
  !!
  subroutine pdftab_AllocTab_(gd, tab, Qmin, Qmax, dlnlnQ, freeze_at_Qmin)
    type(grid_def),    intent(in)  :: gd
    type(pdftab),      intent(out) :: tab
    real(dp), intent(in)           :: Qmin, Qmax
    real(dp), intent(in), optional :: dlnlnQ
    logical,  intent(in), optional :: freeze_at_Qmin
    !----------------------------------------------
    integer :: iQ
    tab%gd = gd
    tab%lambda_eff = min(half*Qmin, default_lambda_eff)
    tab%lnlnQ_min = lnln(tab,Qmin)
    tab%lnlnQ_max = lnln(tab,Qmax)

    tab%default_dlnlnQ = default_or_opt(default_dlnlnQ, dlnlnQ)
    tab%nQ = ceiling((tab%lnlnQ_max - tab%lnlnQ_min)/tab%default_dlnlnQ)
    tab%dlnlnQ = (tab%lnlnQ_max - tab%lnlnQ_min)/tab%nQ
    tab%freeze_at_Qmin = default_or_opt(.false.,freeze_at_Qmin)

    !-- by default, no extra info is given
    tab%nf_info_associated = .false.
    nullify(tab%as2pi)
    nullify(tab%nf_int)
    nullify(tab%evops)

    !write(0,*) 'pdftab info: Number of Q bins is ',tab%nQ
    call pdfgen_AllocPDF(gd,tab%tab,0,tab%nQ)
    allocate(tab%lnlnQ_vals(0:tab%nQ))
    do iQ = 0, tab%nQ
       tab%lnlnQ_vals(iQ) = tab%lnlnQ_min + iQ * tab%dlnlnQ
    end do
    
  end subroutine pdftab_AllocTab_


  
  !---------------------------------------------------------
  !! 1d overloaded version of pdftab_AllocTab
  subroutine pdftab_AllocTab_1d(gd, tab, Qmin, Qmax, dlnlnQ, freeze_at_Qmin)
    type(grid_def), intent(in)     :: gd
    type(pdftab),   intent(out)    :: tab(:)
    real(dp), intent(in)           :: Qmin, Qmax
    real(dp), intent(in), optional :: dlnlnQ
    logical,  intent(in), optional :: freeze_at_Qmin
    integer :: i

    do i = 1, size(tab)
       call pdftab_AllocTab_(gd, tab(i), Qmin, Qmax, dlnlnQ, freeze_at_Qmin)
    end do
  end subroutine pdftab_AllocTab_1d


  !---------------------------------------------------------------
  !! Associate a tab with the nf info from an alpha_s holder
  !! (i.e. read that info, and copy it into tab).
  !!
  !! This involves REALLOCATING the tab pointer so as to allow for
  !! extra information.
  !!
  !! The new tab separates stretches with different nf values. It also
  !! has extra info in the form of arrays nf_int and as2pi, which make
  !! it possible to write expressions such as
  !!
  !!   Dtab%tab = tab%as2pi(:) * (sh%P(iloop,tab%nf_int(:)) .conv. tab%tab)
  !!
  !! This is one of the rare cases in which direct access to structure 
  !! components is still allowed... (Lack of finalisation makes it
  !! difficult to do otherwise).
  !!
  !! KNOWN LIMITATIONS: What happens with muM_mQ /= 1? All hell breaks
  !! loose?
  subroutine pdftab_AssocNfInfo(tab,ash)
    type(pdftab), intent(inout) :: tab
    type(as_handle), intent(in) :: ash
    !-----------------------------------
    integer  :: nflcl, iQ_prev, iQ
    real(dp) :: Qlo, Qhi, Qhi_test
    type(pdfseginfo), pointer :: seginfo

    
    if (tab%nf_info_associated) call wae_error('pdftab_AssocNfInfo',&
         &'nf info already associated: delete it first')

    ! We will be reallocating everything here, so first clean up
    call pdftab_DelTab(tab)
    tab%dlnlnQ = zero ! will no longer be useful...

    call as_nfatQ(ash, invlnln(tab,tab%lnlnQ_min), tab%nflo)
    call as_nfatQ(ash, invlnln(tab,tab%lnlnQ_max), tab%nfhi)

    allocate(tab%seginfo(tab%nflo:tab%nfhi))

    ! figure out how we are going to bin things...
    iQ_prev = -1
    do nflcl = tab%nflo, tab%nfhi
       !write(0,*) 'nflcl',nflcl
       seginfo => tab%seginfo(nflcl)

       call as_QRangeAtNf(ash, nflcl, Qlo, Qhi_test)
       call as_QRangeAtNf(ash, nflcl, Qlo, Qhi, muM_mQ=one)
       ! if one weakens this restriction, then one should think about
       ! the consequences for the determination of alpha_s/2pi, here
       ! and elsewhere....
       if (Qhi_test /= Qhi) call wae_error('pdftab_AssocNfInfo',&
         &'it seems that ash has muM_mQ /= one. Currently unsupported.',&
         &dbleval=Qhi_test/Qhi)

       ! Include min_dlnlnQ_singleQ to ensure we do not go EXACTLY to the
       ! mass threshold, where, in evolution one might run into problems.
       ! BUT, we will later have to worry about what to do when we
       ! are in between thresholds...
       seginfo%lnlnQ_lo = max(lnln(tab,Qlo)+min_dlnlnQ_singleQ, tab%lnlnQ_min)
       seginfo%lnlnQ_hi = min(lnln(tab,Qhi)-min_dlnlnQ_singleQ, tab%lnlnQ_max)

       seginfo%ilnlnQ_lo = iQ_prev + 1
       !write(0,*) 'ill_lo', seginfo%ilnlnQ_lo
       if ((seginfo%lnlnQ_hi - seginfo%lnlnQ_lo) < two*min_dlnlnQ_singleQ) then
          ! use just one point
          seginfo%ilnlnQ_hi = seginfo%ilnlnQ_lo
          seginfo%dlnlnQ = zero
          seginfo%lnlnQ_hi = seginfo%lnlnQ_lo
       else
          seginfo%ilnlnQ_hi = seginfo%ilnlnQ_lo + max(lnlnQ_order,&
               & ceiling((seginfo%lnlnQ_hi - seginfo%lnlnQ_lo)/&
               & tab%default_dlnlnQ))
          seginfo%dlnlnQ = (seginfo%lnlnQ_hi - seginfo%lnlnQ_lo)/&
               & (seginfo%ilnlnQ_hi - seginfo%ilnlnQ_lo)
       end if
       !write(0,*) 'ill_hi', seginfo%ilnlnQ_hi, seginfo%dlnlnQ, &
       !     &invlnln(tab,seginfo%lnlnQ_lo),invlnln(tab,seginfo%lnlnQ_hi), tab%default_dlnlnQ
       iQ_prev = seginfo%ilnlnQ_hi
    end do
    
    ! this should not happen too often! But check it just in
    ! case...
    if (tab%seginfo(tab%nflo)%lnlnQ_lo /= tab%lnlnQ_min .or.&
         &tab%seginfo(tab%nfhi)%lnlnQ_hi /= tab%lnlnQ_max) &
         &call wae_error('pdftab_AssocNfInfo',&
         & 'mismatch in segment and global lnlnQ limits.',&
         & 'Could be due ash having more restricted range?')


    ! now reallocate things?
    tab%nQ = tab%seginfo(tab%nfhi)%ilnlnQ_hi
    call pdfgen_AllocPDF(tab%gd,tab%tab,0,tab%nQ)
    allocate(tab%lnlnQ_vals(0:tab%nQ))
    allocate(tab%nf_int(0:tab%nQ))
    allocate(tab%as2pi(0:tab%nQ))

    ! set up complementary info...
    do nflcl = tab%nflo, tab%nfhi
       !write(0,*) 'nflcl',nflcl
       seginfo => tab%seginfo(nflcl)
       do iQ = seginfo%ilnlnQ_lo, seginfo%ilnlnQ_hi
          tab%nf_int(iQ) = nflcl
          tab%lnlnQ_vals(iQ) = seginfo%lnlnQ_lo &
               & + (iQ-seginfo%ilnlnQ_lo)*seginfo%dlnlnQ
          tab%as2pi(iQ) = as_Value(ash,invlnln(tab,tab%lnlnQ_vals(iQ)))/twopi
       end do
    end do
    
    ! REMEMBER TO COMPLETE FROM ORIG...
    tab%nf_info_associated = .true.
    write(0,*) 'pdftab info: Number of Q bins changed to',tab%nQ
  end subroutine pdftab_AssocNfInfo


  !---------------------------------------------------------------
  !! 1d-overloaded verseion of pdftab_AssocNfInfo
  subroutine pdftab_AssocNfInfo_1d(tab,ash)
    type(pdftab), intent(inout) :: tab(:)
    type(as_handle), intent(in) :: ash
    !-----------------------------------
    integer :: i
    do i = 1, size(tab)
       call pdftab_AssocNfInfo(tab(i),ash)
    end do
    
  end subroutine pdftab_AssocNfInfo_1d


  !-----------------------------------------------------------------------
  !! Allocate the memory for a new tab, using as a template an
  !! preexistent tab (origtab).
  !!
  !! Additionally, information concerning any varnf and alphas
  !! structure copied from origtab to tab. Actual PDF contents of the
  !! tab are not however copied.
  !! 
  subroutine pdftab_AllocTab_fromorig(tab, origtab)
    type(pdftab), intent(out) :: tab
    type(pdftab), intent(in)  :: origtab

    tab = origtab
    !-- this is the only thing that is not taken care of...
    call pdfgen_AllocPDF(tab%gd,tab%tab,0,tab%nQ)
    allocate(tab%lnlnQ_vals(0:tab%nQ))
    tab%lnlnQ_vals = origtab%lnlnQ_vals

    if (origtab%nf_info_associated) then
       allocate(tab%seginfo(tab%nflo:tab%nfhi))
       allocate(tab%nf_int(0:tab%nQ))
       allocate(tab%as2pi(0:tab%nQ))
       tab%seginfo = origtab%seginfo
       tab%nf_int = origtab%nf_int
       tab%as2pi = origtab%as2pi
    end if
  end subroutine pdftab_AllocTab_fromorig

  !---------------------------------------------------------
  !! 1d-overloaded version of pdftab_AllocTab_fromorig
  subroutine pdftab_AllocTab_fromorig_1d(tab, origtab)
    type(pdftab), intent(out) :: tab(:)
    type(pdftab), intent(in)  :: origtab
    integer :: i
    do i = 1, size(tab)
       call pdftab_AllocTab_fromorig(tab(i), origtab)
    end do
  end subroutine pdftab_AllocTab_fromorig_1d


  !------------------------------------------------------
  !! Initialise the tab with the results of a subroutine (see
  !! interface). Note that the subroutine takes y = ln1/x and Q as its
  !! arguments.
  subroutine pdftab_InitTabSub_(tab, sub)
    type(pdftab), intent(inout) :: tab
    interface
       subroutine sub(y, Q, res)
         use types; implicit none
         real(dp), intent(in) :: y, Q
         real(dp), intent(out):: res(:)
       end subroutine sub
    end interface
    !-------------
    integer :: iQ
    real(dp) :: Q
    do iQ = 0, tab%nQ
       Q = invlnln(tab,tab%lnlnQ_min + iQ*tab%dlnlnQ)
       call pdfgen_InitPDFSub(tab%gd, tab%tab(:,:,iQ), sub, Q)
    end do
    
  end subroutine pdftab_InitTabSub_
  
  !------------------------------------------------------
  !! Initialise the tab with the results of a subroutine (see
  !! interface). In addition y = ln1/x and Q, the subroutine takes the
  !! argument iset, enabling one to initialise from a subroutine that
  !! provides several "PDF sets".
  !!
  subroutine pdftab_InitTabSub_iset(tab, sub, iset)
    type(pdftab), intent(inout) :: tab
    integer,      intent(in)    :: iset
    interface
       subroutine sub(y, Q, iset, res)
         use types; implicit none
         real(dp), intent(in) :: y, Q
         integer, intent(in) :: iset
         real(dp), intent(out):: res(:)
       end subroutine sub
    end interface
    !-------------
    integer :: iQ
    real(dp) :: Q
    do iQ = 0, tab%nQ
       !Q = invlnln(tab,tab%lnlnQ_min + iQ*tab%dlnlnQ)
       Q = invlnln(tab,tab%lnlnQ_vals(iQ))
       call pdfgen_InitPDFSub(tab%gd, tab%tab(:,:,iQ), sub, Q, iset)
    end do
    
  end subroutine pdftab_InitTabSub_iset

  !------------------------------------------------------
  !! Initialise the tab with the results of a subroutine (see
  !! interface). Note that the subroutine takes y = ln1/x and Q as its
  !! arguments.
  subroutine pdftab_InitTab_LHAPDF(tab, LHAsub)
    type(pdftab), intent(inout) :: tab
    interface
       subroutine LHAsub(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine LHAsub
    end interface
    !-------------
    integer :: iQ
    real(dp) :: Q
    do iQ = 0, tab%nQ
       Q = invlnln(tab,tab%lnlnQ_min + iQ*tab%dlnlnQ)
       call pdfgen_InitPDF_LHAPDF(tab%gd, tab%tab(:,:,iQ), LHAsub, Q)
    end do
  end subroutine pdftab_InitTab_LHAPDF

  !---------------------------------------------------------------------
  !! Given a starting distribution, StartDist, at StartScale and an
  !! "ash" determining the behaviour of alphas, fill in the tab with
  !! the evolution of the starting distribution. Most of the arguments
  !! have meanings similar to those in the evolution routines
  !!
  subroutine pdftab_InitTabEvolve(tab, StartScale, StartDist, &
                                            & sh, nloop, ash, muR_Q,untie_nf)
    type(pdftab),    intent(inout) :: tab
    real(dp),           intent(in) :: StartScale
    real(dp),           intent(in) :: StartDist(0:,id_min:)
    type(sigma_holder), intent(in) :: sh
    integer,            intent(in) :: nloop
    type(as_handle),    intent(in) :: ash
    real(dp), intent(in), optional :: muR_Q
    logical,  intent(in), optional :: untie_nf
    !-----------------------------------------------------
    
    call pdftab_TabEvolveGen(tab, StartScale, sh, nloop, ash, &
       &StartDist=StartDist, muR_Q=muR_Q, untie_nf=untie_nf)
  end subroutine pdftab_InitTabEvolve
  
  !---------------------------------------------------------------------
  !! Given a starting scale, precalculate the evolution operators that
  !! are needed to generate the full tab from a distribution at that 
  !! starting scale . Most of the arguments have meanings similar
  !! to those in the evolution routines
  !!
  subroutine pdftab_PreEvolve(tab, StartScale, sh, nloop, ash, muR_Q,untie_nf)
    type(pdftab),    intent(inout) :: tab
    real(dp),           intent(in) :: StartScale
    type(sigma_holder), intent(in) :: sh
    integer,            intent(in) :: nloop
    type(as_handle),    intent(in) :: ash
    real(dp), intent(in), optional :: muR_Q
    logical,  intent(in), optional :: untie_nf
    !-----------------------------------------------------
    
    call pdftab_TabEvolveGen(tab, StartScale, sh, nloop, ash, &
       &precalc = .true., muR_Q=muR_Q, untie_nf=untie_nf)
  end subroutine pdftab_PreEvolve

  !---------------------------------------------------------------------
  !! Given a starting distribution, StartDist, and assuming that the
  !! pre-evolution has been carried out for the tab, then generate the
  !! contents of the table from the StartDist.
  !!
  !! NB: alphas will not be updated relative to last true evolution or
  !!     preevolution. So if at any point since doing the preevolution,
  !!     an evolution has been carried out with a different alphas, then 
  !!     as2pi will not be correct here. SHOULD THIS BE FIXED?
  !!
  subroutine pdftab_InitTabEvolve_frompre(tab, StartDist)
    type(pdftab),    intent(inout) :: tab
    real(dp),           intent(in) :: StartDist(0:,id_min:)
    !-----------------------------------------------------
    real(dp) :: dist(0:ubound(StartDist,1),id_min:ubound(StartDist,2))
    integer :: iQ
    
    if (.not. associated(tab%evops)) call wae_error(&
         &'pdftab_InitTabEvolve_frompre',&
         &'No precalculated evolution is available')

    dist = StartDist
    do iQ = tab%StartScale_iQlo, 0, -1
       dist = tab%evops(iQ) .conv. dist
       tab%tab(:,:,iQ) = dist
    end do
    dist = StartDist
    do iQ = tab%StartScale_iQlo+1, tab%nQ
       dist = tab%evops(iQ) .conv. dist
       tab%tab(:,:,iQ) = dist
    end do
    
  end subroutine pdftab_InitTabEvolve_frompre

  !---------------------------------------------------------------------
  !! General internal routine, which serves both to carry out evolution
  !! of a parton distribution AND to determine the evops. That way
  !! everything to do with evolution is concentrated in a single location.
  !!
  !! For the user, this routine should probably be accessed via the more
  !! specific routines, pdftab_InitTabEvolve and pdftab_PreCalc.
  !!
  !! Given a starting distribution, StartDist, at StartScale and an
  !! "ash" determining the behaviour of alphas, fill in the tab with
  !! the evolution of the starting distribution.
  !!
  subroutine pdftab_TabEvolveGen(tab, StartScale, sh, nloop, ash, &
       &StartDist, precalc, muR_Q,untie_nf)
    type(pdftab),    intent(inout) :: tab
    real(dp),           intent(in) :: StartScale
    type(sigma_holder), intent(in) :: sh
    integer,            intent(in) :: nloop
    type(as_handle),    intent(in) :: ash
    !real(dp), optional, intent(in) :: StartDist(0:,id_min:)
    real(dp), optional, intent(in) :: StartDist(:,:)
    logical,  optional, intent(in) :: precalc
    real(dp), intent(in), optional :: muR_Q
    logical,         intent(in), optional  :: untie_nf
    !-----------------------------------------------------
    real(dp), allocatable :: dist(:,:)
    real(dp) :: lnlnQ_norm, lnlnQ, Q_init, Q_end, last_Q
    integer :: i, iQ_lo, iQ_hi
    logical :: precalc_lcl

    precalc_lcl = default_or_opt(.false.,precalc)
    if (precalc_lcl) then
       if (associated(tab%evops)) call wae_error('pdftab_TabEvolveGen',&
            &'tab%evops has already been calculated. Delete the tab first,',&
            &'if you want to recalculated it.')
       allocate(tab%evops(0:tab%nQ))
    end if
    
    
    lnlnQ = lnln(tab,StartScale)
    call request_iQrange(tab, lnlnQ, 1, iQ_lo, iQ_hi, lnlnQ_norm)
    tab%StartScale = StartScale
    tab%StartScale_iQlo = iQ_lo
    ! force this...
    iQ_hi = iQ_lo + 1
    write(0,*) iQ_lo, iQ_hi
    !lnlnQ_norm_start = lnln_norm(tab,StartScale)


    if (present(StartDist)) then
       allocate(dist(size(StartDist,1),size(StartDist,2)))
       dist = StartDist
    end if
    
    last_Q = StartScale
    !do i = floor(lnlnQ_norm_start), 0, -1
    do i = iQ_lo, 0, -1
       Q_init = last_Q
       Q_end = invlnln(tab,tab%lnlnQ_vals(i))
       !write(0,*) 'doing ev from ',Q_init,' to', Q_end

       if (present(StartDist)) then
          call ev_evolve_varnf(sh, dist, &
               &nloop, ash, Q_init, Q_end, muR_Q, untie_nf)
          tab%tab(:,:,i) = dist
       end if
       
       if (precalc_lcl) call ev_evolve_varnf_gen(sh, &
            &nloop, ash, Q_init, Q_end, evop=tab%evops(i), &
            &muR_Q = muR_Q, untie_nf = untie_nf)

       if (tab%nf_info_associated) tab%as2pi(i) = as_Value(ash,Q_end)/twopi
       last_Q = Q_end
    end do
    

    if (present(StartDist)) dist = StartDist
    last_Q = StartScale
    !do i = ceiling(lnlnQ_norm_start), tab%nQ
    do i = iQ_hi, tab%nQ
       Q_init = last_Q
       Q_end = invlnln(tab,tab%lnlnQ_vals(i))
       !write(0,*) 'doing ev from ',Q_init,' to', Q_end

       !write(0,*) 'doing ev from ',Q_init,' to', Q_end
       if (present(StartDist)) then
          call ev_evolve_varnf(sh, dist, nloop,&
               & ash, Q_init, Q_end, muR_Q, untie_nf)
          tab%tab(:,:,i) = dist
       end if
       
       if (precalc_lcl) call ev_evolve_varnf_gen(sh, &
            &nloop, ash, Q_init, Q_end, evop=tab%evops(i), &
            &muR_Q = muR_Q, untie_nf = untie_nf)

       if (tab%nf_info_associated) tab%as2pi(i) = as_Value(ash,Q_end)/twopi
       last_Q = Q_end
    end do
    
    if (present(StartDist)) deallocate(dist)
  end subroutine pdftab_TabEvolveGen


  !--------------------------------------------------------------------
  !! Returns a vector val(id_min:id_max) for the PDF at this
  !! y=ln1/x,Q.
  subroutine pdftab_ValTab_yQ(tab,y,Q,val)
    type(pdftab), intent(in) :: tab
    real(dp),     intent(in) :: y, Q
    real(dp),    intent(out) :: val(id_min:)
    !----------------------------------------
    real(dp) :: lnlnQ, lnlnQ_norm
    real(dp) :: lnlnQ_wgts(0:lnlnQ_order)
    real(dp), pointer :: y_wgts(:)
    real(dp), allocatable :: wgts(:,:)
    integer :: ilnlnQ_lo, ilnlnQ_hi, nQ,iylo, iQ, if
    integer, save :: warn_id = warn_id_INIT

    if (ubound(val,dim=1) < id_max) call wae_error('pdftab_ValTab',&
         &'upper bound of val is too low', intval=ubound(val,dim=1))

    !-- y weights taken care of elsewhere....
    call conv_WgtGridQuant(tab%gd, y, iylo, y_wgts)
    !-- Q weights need some help in finding location etc.
    lnlnQ = lnln(tab,Q)

    if (tab%freeze_at_Qmin .and. lnlnQ < tab%lnlnQ_min) then
       lnlnQ = tab%lnlnQ_min
    else if (lnlnQ < (one-warn_tolerance)*tab%lnlnQ_min &
         &.or. lnlnQ > (one+warn_tolerance)*tab%lnlnQ_max) then
       call wae_warn(default_max_warn, warn_id, &
            &'pdftab_ValTab: Q out of range; result set to zero; Q was:',&
            &dbleval=Q)
       val = zero
       return
    end if
    
    call request_iQrange(tab,lnlnQ,lnlnQ_order,ilnlnQ_lo,ilnlnQ_hi,lnlnQ_norm)
    nQ = ilnlnQ_hi - ilnlnQ_lo
    !OLD lnlnQ_norm = (lnlnQ-tab%lnlnQ_min)/tab%dlnlnQ
    !OLD ilnlnQ_lo = floor(lnlnQ_norm) - lnlnQ_order/2
    !OLD ilnlnQ_lo = max(0, min(tab%nQ-lnlnQ_order,ilnlnQ_lo))
    call uniform_interpolation_weights(lnlnQ_norm, lnlnQ_wgts(0:nQ))

    allocate(wgts(lbound(y_wgts,dim=1):ubound(y_wgts,dim=1),0:nQ))
    do iQ = 0, nQ
       wgts(:,iQ) = y_wgts(:) * lnlnQ_wgts(iQ)
    end do

    !-- is this order more efficient, or should we not bother to
    !   calculate wgts?
    do if = id_min, id_max
       val(if) = sum(wgts*tab%tab(iylo:iylo+size(y_wgts)-1,&
            & if,ilnlnQ_lo:ilnlnQ_hi))
    end do
    !write(0,*) ilnlnQ_lo, ilnlnQ_hi, real(lnlnQ_wgts), val(1)
    
    deallocate(y_wgts, wgts)
  end subroutine pdftab_ValTab_yQ


  !----------------------------------------------------------------
  !! Returns a vector val(id_min:id_max) for the PDF at this x,Q.
  subroutine pdftab_ValTab_xQ(tab,x,Q,val)
    type(pdftab), intent(in) :: tab
    real(dp),     intent(in) :: x, Q
    real(dp),    intent(out) :: val(id_min:)
    call pdftab_ValTab_yQ(tab,-log(x),Q,val)
  end subroutine pdftab_ValTab_xQ

  
  !-----------------------------------------------------------
  !! Deletes all allocated info associated with the tabulation
  subroutine pdftab_DelTab_0d(tab)
    type(pdftab), intent(inout) :: tab
    integer :: i
    deallocate(tab%tab)
    deallocate(tab%lnlnQ_vals)
    if (tab%nf_info_associated) then
       deallocate(tab%seginfo)
       deallocate(tab%nf_int)
       deallocate(tab%as2pi)
    end if
    if (associated(tab%evops)) then
       do i = 0, tab%nQ
          call ev_DelEvOp(tab%evops(i))
       end do
       !write(0,*) "CURRENTLY UNABLE TO DELETE EVOP CONTENTS"
       deallocate(tab%evops)
    end if
  end subroutine pdftab_DelTab_0d
  subroutine pdftab_DelTab_1d(tab)
    type(pdftab), intent(inout) :: tab(:)
    integer :: i
    do i = 1, size(tab)
       call pdftab_DelTab_0d(tab(i))
    end do
  end subroutine pdftab_DelTab_1d
  

  !------------------------------------------------------------------
  !! Given a tab and lnlnQ value, determine the range of iQ entries,
  !! iQ_lo:iQ_hi to be used for interpolating the grid.
  !!
  !! Where possible (i.e. if there are sufficient Q values in the
  !! grid) ensure that iQ_hi-iQ_lo = nrequest; for grids with varnf,
  !! all iQ_lo:iQ_hi values should correspond to the same nf value.
  !! 
  !! return value of lnlnQ_norm = (lnlnQ - lnlnQ(iQ_lo))/dlnlnQ
  !!
  subroutine request_iQrange(tab, lnlnQ, nrequest, iQ_lo, iQ_hi, lnlnQ_norm)
    type(pdftab), intent(in) :: tab
    real(dp),     intent(in) :: lnlnQ
    integer,      intent(in) :: nrequest
    integer,     intent(out) :: iQ_lo, iQ_hi
    real(dp),    intent(out) :: lnlnQ_norm
    !-------------------------------------
    real(dp) :: dist, distclosest
    integer  :: nfclosest, nflcl
    type(pdfseginfo), pointer :: seginfo
    if (nrequest < 1) call wae_error('request_iQrange',&
         &'nrequest should be >= 1',intval=nrequest)

    if (.not. tab%nf_info_associated) then
       lnlnQ_norm = (lnlnQ-tab%lnlnQ_min)/tab%dlnlnQ
       iQ_lo = floor(lnlnQ_norm) - nrequest/2
       iQ_lo = max(0, min(tab%nQ-nrequest,iQ_lo))
       iQ_hi = min(tab%nQ,iQ_lo + nrequest)
       lnlnQ_norm = lnlnQ_norm - iQ_lo
    else
       ! need to find range to which we are closest. There is certainly
       ! a better way of doing it...
       distclosest = 1e10_dp
       do nflcl = tab%nflo, tab%nfhi
          dist = max(zero,(lnlnQ-tab%seginfo(nflcl)%lnlnQ_hi),&
               & (tab%seginfo(nflcl)%lnlnQ_lo-lnlnQ))
          if (dist < distclosest) then
             nfclosest = nflcl
             distclosest = dist
          end if
       end do
       seginfo => tab%seginfo(nfclosest)
       if (seginfo%ilnlnQ_lo == seginfo%ilnlnQ_hi) then
          iQ_lo = seginfo%ilnlnQ_lo
          iQ_hi = iQ_lo
          lnlnQ_norm = zero
       else
          lnlnQ_norm = (lnlnQ-seginfo%lnlnQ_lo)/seginfo%dlnlnQ &
               &+ seginfo%ilnlnQ_lo
          !write(0,*) lnlnQ_norm, seginfo%dlnlnQ
          iQ_lo = floor(lnlnQ_norm) - nrequest/2
          iQ_lo = max(seginfo%ilnlnQ_lo, min(seginfo%ilnlnQ_hi-nrequest,iQ_lo))
          iQ_hi = min(seginfo%ilnlnQ_hi,iQ_lo + nrequest)
          lnlnQ_norm = lnlnQ_norm - iQ_lo
          !write(0,*) iQ_lo, iQ_hi, invlnln(tab,lnlnQ), lnlnQ_norm
       end if
    end if
  end subroutine request_iQrange
  

  !-------------------------------------
  !! conversion from Q to lnlnQ
  function lnln(tab,Q)
    type(pdftab), intent(in) :: tab
    real(dp),     intent(in) :: Q
    real(dp)                 :: lnln

    lnln = log(log(Q/tab%lambda_eff))
  end function lnln

  !! conversion from lnlnQ to Q
  function invlnln(tab,lnlnQ)
    type(pdftab), intent(in) :: tab
    real(dp),     intent(in) :: lnlnQ
    real(dp)                 :: invlnln

    invlnln = exp(exp(lnlnQ))*tab%lambda_eff
  end function invlnln

!OLD! function lnln_norm(tab,Q)
!OLD! type(pdftab), intent(in) :: tab
!OLD! real(dp), intent(in) :: Q
!OLD! real(dp) :: lnln_norm
!OLD!
!OLD! if (tab%nf_info_associated) call wae_error('lnln_norm')
!OLD! lnln_norm = (lnln(tab,Q) - tab%lnlnQ_min)/tab%dlnlnQ
!OLD! end function lnln_norm
!OLD!
!OLD! function invlnln_norm(tab,lnlnQ_norm)
!OLD! type(pdftab), intent(in) :: tab
!OLD! real(dp), intent(in) :: lnlnQ_norm
!OLD! real(dp) :: invlnln_norm
!OLD!
!OLD! if (tab%nf_info_associated) call wae_error('invlnln_norm')
!OLD! invlnln_norm = invlnln(tab, lnlnQ_norm*tab%dlnlnQ + tab%lnlnQ_min)
!OLD! end function invlnln_norm
end module pdf_tabulate_new
