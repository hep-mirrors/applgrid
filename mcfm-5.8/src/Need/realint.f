      double precision function realint(vector,wgt)
      implicit none
      include 'constants.f'
      include 'debug.f'
      include 'realonly.f'
      include 'virtonly.f'
      include 'noglue.f'
      include 'nflav.f'
      include 'vegas_common.f'
      include 'ptilde.f'
      include 'new_pspace.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'realwt.f'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'process.f'
      include 'PDFerrors.f'
      include 'masses.f'
      include 'wts_bypart.f'
      include 'dipolescale.f'
      include 'stopscales.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'

c P.S.
      include 'APPLinclude.f'
      logical debug222
      integer iProcess
      double precision realContrib(0:maxd,0:nSubProcess-1)
      double precision totalFactor
      double precision fillContrib(0:nSubProcess-1)
      double precision procReference(0:maxd,0:nSubProcess-1)
      double precision currentWeight
c  P.S. end

cz
      integer nproc
      common/nproc/nproc
cz //
      integer ih1,ih2,j,k,nd,nmax,nmin,nvec
      double precision vector(mxdim),W,val,val2,valsum,xint,dot,tmp
      double precision sqrts,fx1(-nf:nf),fx2(-nf:nf),
     . dipfx1(0:maxd,-nf:nf),dipfx2(0:maxd,-nf:nf),
     . fx1_H(-nf:nf),fx2_H(-nf:nf),fx1_L(-nf:nf),fx2_L(-nf:nf)
      double precision p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
      double precision pswt,rscalestart,fscalestart
      double precision s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      double precision msqc(maxd,-nf:nf,-nf:nf),xmsq(0:maxd)
      double precision flux,BrnRat,xreal,xreal2
      double precision xx1,xx2,q(mxpart,4),msqa(-nf:nf,-nf:nf)
      double precision m3,m4,m5,R,Rbbmin
      double precision xmsq_bypart(0:maxd,-1:1,-1:1),xmsqjk,
     . lord_bypart(-1:1,-1:1)
      integer n2,n3,sgnj,sgnk
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
      common/xreal/xreal,xreal2
      common/Rbbmin/Rbbmin
      logical bin,first,failed
      logical incldip(0:maxd),includedipole,includereal
      logical creatent,dswhisto
      character*30 runstring
      common/runstring/runstring
      external qqb_w2jet_g,qqb_w2jet_gs,qqb_z2jet_g,qqb_z2jet_gs,
     . qqb_w2jet,qqb_w1jet_gs,qqb_z2jet,qqb_z1jet_gs,qqb_Hg_g,qqb_Hg_gs,
     . qqb_hww_g,qqb_hww_gs,qqb_zbb_g,qqb_zbb_gs,
     . qqb_wbb_g,qqb_wbb_gs,
     . qqb_w_g,qqb_w_gs,qqb_z1jet,qqb_z_gs,qqb_ww_g,qqb_ww_gs,
     . qqb_wz_g,qqb_wz_gs,qqb_zz_g,qqb_zz_gs,qqb_wgam_g,qqb_wgam_gs,
     . qqb_QQb_g,qqb_QQb_gs,
     . VV_Hqq_g,VV_Hqq_gs,VV_HWW_g,VV_HWW_gs,
     . gg_Hg,gg_H_gs,
     . gg_HWWgg,gg_HWWg_gs,gg_HZZgg,gg_HZZg_gs,
     . gg_Hgg,gg_Hg_gs,
     . gQ_zQ_g,gQ_zQ_gs,qqb_tbb_g,qqb_tbb_gs,
     . qqb_w_tndk_g,qqb_w_tndk_gs,
     . qqb_w_twdk_g,qqb_w_twdk_gs,qqb_w_twdk_gdk,qqb_w_twdk_gsdk,
     . qqb_zbjet_g,qqb_zbjet_gs,qqb_w_cjet_g,qqb_w_cjet_gs,
     . gg_hggg,gg_hgg_gs,
     . qg_tbq_g,qg_tbq_gs,qq_tbg_g,qq_tbg_gs,epem3j_g,epem3j_gs,
     . qqb_QQbdk_g,qqb_QQbdk_gs
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/Pext/p1ext,p2ext
      common/nmax/nmax
      common/BrnRat/BrnRat
      common/nmin/nmin
      common/bypart/lord_bypart
      common/incldip/incldip
      common/outputflags/creatent,dswhisto
cz Add b fraction
      double precision bwgt
      common/btagging/ bwgt
      double precision msqtmp(0:maxd),bwgttmp(0:maxd)
      double precision realeventp(mxpart,4)
      common/realeventp/realeventp
      double precision y32,y43,z3,z4,z5,z6
      
      data bwgt / 0d0 /  ! in common block
cz // Add b fraction   Note: only msqtmp(0), bwgttmp(0) are used in nplotter.f
      data p/48*0d0/
      data first/.true./
      save first,rscalestart,fscalestart
      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif
      ntotshot=ntotshot+1
      pswt=0d0
      realint=0d0      

      W=sqrts**2
      
      if (first) then
         write(6,*)
         write(6,*) 'nmin=',nmin,',nmax=',nmax
         write(6,*)
         first=.false.
      endif
      
c--- processes that use "gen3"     
      if     ( (case .eq. 'W_only')
     .    .or. (case .eq. 'Z_only')
     .    .or. (case .eq. 'ggfus0')
     .   .or.  (case .eq. 'Wcsbar')
     .   .or.  (case .eq. 'Wcs_ms') ) then
        npart=3
        if (new_pspace) then
          call gen3a(vector,p,pswt,*999)
        else
          call gen3(vector,p,pswt,*999)
        endif

c--- processes that use "gen3m"     
      elseif ( (case .eq. 'tt_tot')
     .    .or. (case .eq. 'bb_tot')
     .    .or. (case .eq. 'cc_tot') ) then
        m3=mass2
        m4=mass2
        m5=0d0
        npart=3
        call gen3m(vector,p,m3,m4,m5,pswt,*999)
          
c--- processes that use "gen4"     
      elseif ( (case .eq. 'W_cjet')
     .   .or.  (case .eq. 'Wgamma')
     .   .or.  (case .eq. 'Zgamma')
     .   .or.  (case .eq. 'W_tndk')
     .   .or.  (case .eq. 'epem3j')   ) then
        npart=4
        call gen4(vector,p,pswt,*999)
                  
c--- processes that use "gen4"     
      elseif ((case .eq. 'qg_tbq') .or. (case .eq. 'qq_tbg')) then
        npart=4
        call gen4(vector,p,pswt,*999)
                  
c--- processes that use "gen6"     
      elseif ( 
     .      (case .eq. 'W_twdk') 
     . .or. (case .eq. 'Wtdkay')
     . .or. (case .eq. 'HWWjet')
     . .or. (case .eq. 'HZZjet')
     . ) then
        npart=6
        call gen6(vector,p,pswt,*999)
                  
c--- processes that use "gen7"     
      elseif ( 
     .      (case .eq. 'qq_HWW')
     . .or. (case .eq. 'WH__WW')
     . .or. (case .eq. 'ZH__WW')
     . .or. (case .eq. 'HWW2jt')
     . .or. (case .eq. 'HZZ2jt')
     . ) then
        npart=7
        call gen7(vector,p,pswt,*999)
                  
c--- processes that use "gen7m"     
      elseif ( (case .eq. 'tt_bbl') 
     .    .or. (case .eq. 'tt_bbh')
     . ) then
        m3=mt
        m4=mt
        m5=0d0
        npart=7
	call gen7m(vector,p,m3,m4,m5,pswt,*999)

c--- processes that use "gen_njets" with an argument of "2"     
      elseif ( (case .eq. 'W_1jet')
     .    .or. (case .eq. 'Wcjet0')
     .    .or. (case .eq. 'Z_1jet')
     .    .or. (case .eq. 'H_1jet')
     .    .or. (case .eq. 'ggfus1')
     .    .or. (case .eq. 'gQ__ZQ') ) then
        npart=4
        if (new_pspace) then
          call gen4a(vector,p,pswt,*999)      
        else
          call gen_njets(vector,2,p,pswt,*999)
        endif 
        
c--- processes that use "gen_njets" with an argument of "3"
      elseif ( (case .eq. 'Wbbbar')
     .    .or. (case .eq. 'W_2jet')
     .    .or. (case .eq. 'Z_2jet')
     .    .or. (case .eq. 'Zbbbar')
     .    .or. (case .eq. 'W_bjet')
     .    .or. (case .eq. 'Z_bjet')
     .    .or. (case .eq. 'qq_Hqq')
     .    .or. (case .eq. 'ggfus2') ) then
        npart=5
        call gen_njets(vector,3,p,pswt,*999) 
        
c--- processes that use "gen_stop" with an argument of "1"
      elseif ( (case .eq. 'ttdkay')
     .    .or. (case .eq. 'tdecay') ) then
        npart=5
        call gen_stop(vector,1,p,pswt,*999)

c--- processes that use "gen_stop" with an argument of "2"
      elseif ( (case .eq. 'bq_tpq')
     .    .or. (case .eq. 't_bbar') ) then
        npart=5
        call gen_stop(vector,2,p,pswt,*999)
        
c--- DEFAULT: processes that use "gen5"
      else
        npart=5
        if (new_pspace) then
          call gen5a(vector,p,pswt,*999)
        else
          call gen5(vector,p,pswt,*999)    
        endif
      endif
            
      nvec=npart+2
      call dotem(nvec,p,s)
      
c---impose cuts on final state
      call masscuts(p,*999)
c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
     
c--- extra cut to divide WQj/ZQj regions
      if ( (nproc .eq. 312) .or. (nproc .eq. 317)
     . .or.(nproc .eq. 322) .or. (nproc .eq. 327)
     . .or.(nproc .eq. 342) .or. (nproc .eq. 352)) then
        if (R(p,5,6) .lt. Rbbmin) goto 999
      endif

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead set to zero
      includereal=includedipole(0,p)
      incldip(0)=includereal 
 
      if (includereal .eqv. .false.) then
        do j=-nf,nf
        do k=-nf,nf
          msq(j,k)=0d0
          msqLH(j,k)=0d0   ! for stop+b process
          msqHL(j,k)=0d0   ! for stop+b process
        enddo
        enddo
      endif
      
      do j=1,mxpart
      do k=1,4
        realeventp(j,k)=p(j,k)
      enddo
      enddo
      
      if (dynamicscale) then
        call scaleset(rscalestart,fscalestart,p)
	dipscale(0)=facscale
      endif
      
c---- generate collinear points that satisfy the jet cuts (for checking)
c      call singgen(p,s,*998)
            
c----calculate the x's for the incoming partons from generated momenta

      xx1=two*(p(1,4)*p2ext(4)-p(1,3)*p2ext(3))/W
      xx2=two*(p(2,4)*p1ext(4)-p(2,3)*p1ext(3))/W

      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx1,xx2
      
      if ((xx1 .gt. 1d0) .or. (xx2 .gt. 1d0)) then
         realint=0d0
         return
      endif

c--- Calculate the required matrix elements      
      if     (case .eq. 'W_only') then
c        call singcheck(qqb_w_g,qqb_w_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_w_g(p,msq)      
        call qqb_w_gs(p,msqc)     
      elseif (case .eq. 'W_1jet') then
c        call singcheck(qqb_w2jet,qqb_w1jet_gs,p)   ! Checked 11/16/01
        if (includereal) call qqb_w2jet(p,msq)      
        call qqb_w1jet_gs(p,msqc)  
      elseif (case .eq. 'Wgamma') then
c        call singcheck(qqb_wgam_g,qqb_wgam_gs,p)   ! Checked 08/27/02
        if (includereal) call qqb_wgam_g(p,msq)      
        call qqb_wgam_gs(p,msqc)  
      elseif (case .eq. 'W_cjet') then
c        call singcheck(qqb_w_cjet_g,qqb_w_cjet_gs,p) ! Checked 15/05/07
        if (includereal) call qqb_w_cjet_g(p,msq)      
        call qqb_w_cjet_gs(p,msqc)  
      elseif (case .eq. 'Zgamma') then
c        call singcheck(qqb_zgam_g,qqb_zgam_gs,p)
        if (includereal) call qqb_zgam_g(p,msq)      
        call qqb_zgam_gs(p,msqc)  
      elseif (case .eq. 'Wbbbar') then
c        call singcheck(qqb_wbb_g,qqb_wbb_gs,p)     ! Checked 11/30/01
        if (includereal) call qqb_wbb_g(p,msq)      
        call qqb_wbb_gs(p,msqc)      
      elseif (case .eq. 'W_2jet') then      
c        call singcheck(qqb_w2jet_g,qqb_w2jet_gs,p) ! Re-checked June 09
        if (includereal)  call qqb_w2jet_g(p,msq)  
        call qqb_w2jet_gs(p,msqc)
      elseif (case .eq. 'Z_only') then
c        call singcheck(qqb_z1jet,qqb_z_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_z1jet(p,msq)      
        call qqb_z_gs(p,msqc)     
      elseif (case .eq. 'Z_1jet') then
c        call singcheck(qqb_z2jet,qqb_z1jet_gs,p)   ! Checked 11/16/01
        if (includereal) call qqb_z2jet(p,msq)      
        call qqb_z1jet_gs(p,msqc)  
      elseif (case .eq. 'Z_2jet') then
c        call singcheck(qqb_z2jet_g,qqb_z2jet_gs,p) ! Re-checked June 09
        if (includereal) call qqb_z2jet_g(p,msq)  
        call qqb_z2jet_gs(p,msqc) 
      elseif (case .eq. 'Zbbbar') then
c        call singcheck(qqb_zbb_g,qqb_zbb_gs,p)     ! Checked 11/30/01
        if (includereal) call qqb_zbb_g(p,msq)
        call qqb_zbb_gs(p,msqc) 
      elseif (case .eq. 'WWqqbr') then
c        call singcheck(qqb_ww_g,qqb_ww_gs,p)       ! Checked 11/30/01
        if (includereal) call qqb_ww_g(p,msq)      
        call qqb_ww_gs(p,msqc)      
      elseif (case .eq. 'WZbbar') then
c        call singcheck(qqb_wz_g,qqb_wz_gs,p)       ! Checked 12/05/01
        if (includereal) call qqb_wz_g(p,msq)      
        call qqb_wz_gs(p,msqc)      
      elseif (case .eq. 'ZZlept') then
c        call singcheck(qqb_zz_g,qqb_zz_gs,p)       ! Checked 12/05/01
        if (includereal) call qqb_zz_g(p,msq)      
        call qqb_zz_gs(p,msqc)      
      elseif (case .eq. 'WHbbar') then
c        call singcheck(qqb_wh_g,qqb_wh_gs,p)
        if (includereal) call qqb_wh_g(p,msq)      
        call qqb_wh_gs(p,msqc)     
      elseif (case .eq. 'WH__WW') then
c        call singcheck(qqb_wh_ww_g,qqb_wh_ww_gs,p)
        if (includereal) call qqb_wh_ww_g(p,msq)      
        call qqb_wh_ww_gs(p,msqc)     
      elseif (case .eq. 'ZHbbar') then
c        call singcheck(qqb_zh_g,qqb_zh_gs,p)
        if (includereal) call qqb_zh_g(p,msq)      
        call qqb_zh_gs(p,msqc)     
      elseif (case .eq. 'ZH__WW') then
c        call singcheck(qqb_zh_ww_g,qqb_zh_ww_gs,p)
        if (includereal) call qqb_zh_ww_g(p,msq)      
        call qqb_zh_ww_gs(p,msqc)     
       elseif (case .eq. 'ggfus0') then
c         call singcheck(gg_hg,gg_h_gs,p)       ! Checked 28/02/03
         if (includereal) call gg_hg(p,msq)
         call gg_h_gs(p,msqc)
      elseif (case .eq. 'HWW_4l') then
c        call singcheck(qqb_hww_g,qqb_hww_gs,p)
        if (includereal) call qqb_hww_g(p,msq)      
        call qqb_hww_gs(p,msqc)     
      elseif (case .eq. 'HZZ_4l') then
c        call singcheck(qqb_hzz_g,qqb_hzz_gs,p)
        if (includereal) call qqb_hzz_g(p,msq)      
        call qqb_hzz_gs(p,msqc)  
      elseif (case .eq. 'H_1jet') then
c        call singcheck(qqb_Hg_g,qqb_Hg_gs,p)       ! Checked 19/02/02
        if (includereal) call qqb_Hg_g(p,msq)  
        call qqb_Hg_gs(p,msqc) 
      elseif ((case .eq. 'tt_bbl') .or. (case .eq. 'tt_bbh')) then
c         call singcheck(qqb_QQbdk_g,qqb_QQbdk_gs,p) ! Checked 15/8/08
        if (includereal) call qqb_QQbdk_g(p,msq)  
        call qqb_QQbdk_gs(p,msqc) 
      elseif ((case .eq. 'tt_tot') .or. (case .eq. 'cc_tot')
     .   .or. (case .eq. 'bb_tot')) then
c        call singcheck(qqb_QQb_g,qqb_QQb_gs,p)
        if (includereal) call qqb_QQb_g(p,msq)
        call qqb_QQb_gs(p,msqc)
      elseif ((case .eq. 'bq_tpq') .or. (case .eq. 't_bbar')) then
c        call singcheck(qqb_tbb_g,qqb_tbb_gs,p)
        if (includereal) call qqb_tbb_g(p,msq)
        call qqb_tbb_gs(p,msqc)
      elseif (case .eq. 'ttdkay') then
c        call singcheck(qqb_tbb_g_new,qqb_tbb_gs,p)
        if (includereal) call bq_tpq_gdk(p,msq)
        call bq_tpq_gsdk(p,msqc)
      elseif (case .eq. 'tdecay') then
c        call singcheck(qqb_tbb_g_new,qqb_tbb_gs,p)
        if (includereal)  call qqb_tbb_gdk(p,msq)
        call qqb_tbb_gsdk(p,msqc)       
       elseif (case .eq. 'W_tndk') then
c        call singcheck(qqb_w_tndk_g,qqb_w_tndk_gs,p)	! Checked 12/3/04
        if (includereal) call qqb_w_tndk_g(p,msq)
        call qqb_w_tndk_gs(p,msqc)
      elseif (case .eq. 'W_twdk') then
c        call singcheck(qqb_w_twdk_g,qqb_w_twdk_gs,p)		! Checked 2/4/05
        if (includereal) call qqb_w_twdk_g(p,msq)
        call qqb_w_twdk_gs(p,msqc)
      elseif (case .eq. 'Wtdkay') then
c        call singcheck(qqb_w_twdk_gdk,qqb_w_twdk_gsdk,p)	! Checked 2/2/05
        if (includereal) call qqb_w_twdk_gdk(p,msq)
        call qqb_w_twdk_gsdk(p,msqc)
      elseif (case .eq. 'ggfus1') then
c        call singcheck(gg_hgg,gg_hg_gs,p)
        if (includereal) call gg_hgg(p,msq)
        call gg_hg_gs(p,msqc)
      elseif (case .eq. 'HWWjet') then
c        call singcheck(gg_hWWgg,gg_hWWg_gs,p)
        if (includereal) call gg_hWWgg(p,msq)
        call gg_hWWg_gs(p,msqc)
      elseif (case .eq. 'HWW2jt') then
c        call singcheck(gg_hWWgg,gg_hWWg_gs,p)
        if (includereal) call gg_hWWggg(p,msq)
        call gg_hWWgg_gs(p,msqc)
      elseif (case .eq. 'HZZjet') then
c        call singcheck(gg_hZZgg,gg_hZZg_gs,p)
        if (includereal) call gg_hZZgg(p,msq)
        call gg_hZZg_gs(p,msqc)
      elseif (case .eq. 'HZZ2jt') then
c        call singcheck(gg_hZZgg,gg_hZZg_gs,p)
        if (includereal) call gg_hZZggg(p,msq)
        call gg_hZZgg_gs(p,msqc)
c	write(6,*) msq
c	pause
      elseif (case .eq. 'qq_Hqq') then
c        call singcheck(VV_Hqq_g,VV_Hqq_gs,p)
        if (includereal) call VV_Hqq_g(p,msq)
        call VV_Hqq_gs(p,msqc)
      elseif (case .eq. 'qq_HWW') then
c        call singcheck(VV_HWW_g,VV_HWW_gs,p)
        if (includereal) call VV_HWW_g(p,msq)
        call VV_HWW_gs(p,msqc)
      elseif (case .eq. 'ggfus2') then
c        call singcheck(gg_hggg,gg_hgg_gs,p)		! Checked 10/29/09
        if (includereal) call gg_hggg(p,msq)
        call gg_hgg_gs(p,msqc)
      elseif (case .eq. 'qg_tbq') then
c        call singcheck(qg_tbq_g,qg_tbq_gs,p)	 ! Checked 2/4/08
       if (includereal) call qg_tbq_g(p,msq)
       call qg_tbq_gs(p,msqc)
      elseif (case .eq. 'qq_tbg') then
c        call singcheck(qq_tbg_g,qq_tbg_gs,p)	 ! Checked 8/9/08
       if (includereal) call qq_tbg_g(p,msq)
       call qq_tbg_gs(p,msqc)
      elseif (case .eq. 'epem3j') then
c        call singcheck(epem3j_g,epem3j_gs,p)	 ! Checked 17/11/08
       if (includereal) call epem3j_g(p,msq)
       call epem3j_gs(p,msqc)
      elseif (case .eq. 'gQ__ZQ') then
c        call singcheck(gQ_zQ_g,gQ_zQ_gs,p)
        if (includereal) call gQ_zQ_g(p,msq)
        call gQ_zQ_gs(p,msqc)
      elseif (case .eq. 'Z_bjet') then
c        call singcheck(qqb_zbjet_g,qqb_zbjet_gs,p)	! Checked 07/18/05
        if (includereal) call qqb_zbjet_g(p,msq)
        call qqb_zbjet_gs(p,msqc)
      elseif (case .eq. 'Wcsbar') then
c        call singcheck(qqb_w_g,qqb_w_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_w_g(p,msq)      
        call qqb_w_gs(p,msqc)     
      elseif (case .eq. 'Wcs_ms') then
        if (includereal) call qqb_w_cjet(p,msq)
        ndmax=0
      endif
      
      do nd=0,ndmax
      xmsq(nd)=0d0
cz
      msqtmp(nd)=0d0
      bwgttmp(nd)=0d0
cz //
      do j=-1,1
      do k=-1,1
      xmsq_bypart(nd,j,k)=0d0
      enddo
      enddo
c--->>>P.S. set to zero
      if(creategrid.and.bin)then
         call iniRContribution(nd,realContrib)
         call iniRContribution(nd,procReference)
         debug222 = .false.
      endif
c--->>>P.S. end
      enddo
      
      currentPDF=0
            
      flux=fbGeV2/(two*xx1*xx2*W)
c--- for mlm study, divide by (Ecm)**2=W
      if (runstring(1:3) .eq. 'mlm') then
	flux=flux/W
      endif

c--- initialize a PDF set here, if calculating errors
  777 continue    
      do nd=0,ndmax
      xmsq(nd)=0d0
cz
      msqtmp(nd)=0d0
      bwgttmp(nd)=0d0
cz //
      enddo
      if (PDFerrors) then
        call InitPDF(currentPDF)
      endif
         
c--- calculate PDF's  
      if (dynamicscale) then
        do nd=ndmax,0,-1  ! so that fx1,fx2 correct for real kinematics
          if (dipscale(nd) .lt. 1d-8) then	  
c--- in case dipole is not used, set up dummy value of scale for safety
c--- and set all PDF entries to zero
	    dipscale(nd)=dipscale(0)
	    do j=-nf,nf
	      fx1(j)=0d0
	      fx2(j)=0d0
	    enddo
	  else
            call fdist(ih1,xx1,dipscale(nd),fx1)
            call fdist(ih2,xx2,dipscale(nd),fx2)
            do j=-nf,nf
	      dipfx1(nd,j)=fx1(j)
	      dipfx2(nd,j)=fx2(j)
	    enddo
	  endif
	enddo
      else
        if (case .eq. 'qg_tbq') then
c--- for single top + b, make sure to use two different scales
          call fdist(ih1,xx1,facscale_H,fx1_H)
          call fdist(ih2,xx2,facscale_H,fx2_H)
          call fdist(ih1,xx1,facscale_L,fx1_L)
          call fdist(ih2,xx2,facscale_L,fx2_L)
	  do j=-nf,nf
	    if (j .eq. 0) then  ! heavy quark line has gluon init. state
	      fx1(j)=fx1_H(j)
	      fx2(j)=fx2_H(j)
	    else
	      fx1(j)=fx1_L(j)
	      fx2(j)=fx2_L(j)
	    endif
	  enddo
        else
c--- for comparison with C. Oleari's e+e- --> QQbg calculation
            if (runstring(1:5) .eq. 'carlo') then
            flux=1d0/2d0/W/(as/twopi)**2
c--- divide out by (ason2pi) and then the "LO" massless DY process
	    flux=flux/(aveqq*xn*fourpi*(gwsq/fourpi)**2/3d0/sqrts**2)
  	    flux=flux/(xn/8d0)
	    do j=-nf,nf
	    fx1(j)=0d0
	    fx2(j)=0d0
	    enddo
	    fx1(0)=1d0
	    fx1(1)=1d0
	    fx2(0)=1d0
	    fx2(1)=1d0
          else   
c--- usual case            
            call fdist(ih1,xx1,facscale,fx1)
            call fdist(ih2,xx2,facscale,fx2)
	  endif
        endif
      endif	
            
      do j=-nflav,nflav
      do k=-nflav,nflav

      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (gqonly) then
      if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
      endif      
      
      if (noglue) then 
      if ((j.eq.0) .or. (k.eq.0)) goto 20
      endif

      if ((case .eq. 'Wcsbar').and.(j .ne. 4).and.(k .ne. 4)) goto 20

      if (realonly) then 
        xmsq(0)=xmsq(0)+fx1(j)*fx2(k)*msq(j,k)
        do nd=1,ndmax
        xmsq(nd)=0d0
        enddo
      elseif (virtonly) then
         xmsq(0)=0d0
         do nd=1,ndmax
	   if (dynamicscale) then	   
             xmsq(nd)=xmsq(nd)+dipfx1(nd,j)*dipfx2(nd,k)*(-msqc(nd,j,k))
	   else
             xmsq(nd)=xmsq(nd)+fx1(j)*fx2(k)*(-msqc(nd,j,k))
	   endif
         enddo
      else

         if     (j .gt. 0) then
           sgnj=+1
         elseif (j .lt. 0) then
           sgnj=-1
         else
           sgnj=0
         endif
         if     (k .gt. 0) then
           sgnk=+1
         elseif (k .lt. 0) then
           sgnk=-1
         else
           sgnk=0
         endif

c--- for single top + b, make sure to use two different scales
         if (case .eq. 'qg_tbq') then
           xmsqjk=fx1_L(j)*fx2_H(k)*msqLH(j,k)
     .           +fx1_H(j)*fx2_L(k)*msqHL(j,k)
         else
c--- usual case            
           xmsqjk=fx1(j)*fx2(k)*msq(j,k)
	 endif
	 
         xmsq(0)=xmsq(0)+xmsqjk
cz
cz Extract fraction with b in final state, store in common
cz
c     nproc=161 then t-channel: t
c     isub = 1, nwz = 1
         if (nproc.eq.161) then ! t
            msqtmp(0)=msqtmp(0)+xmsqjk
            if ((j.eq.0).and.((k.lt.0).or.((k.gt.0).and.
     &           (k.ne.5)))) then
               bwgttmp(0)=bwgttmp(0)+xmsqjk
            endif
            if ((k.eq.0).and.((j.lt.0).or.((j.gt.0).and.
     &           (j.ne.5)))) then
               bwgttmp(0)=bwgttmp(0)+xmsqjk
            endif
         endif
c     nproc=166 then t-channel: t~
         if (nproc.eq.166) then ! t
            msqtmp(0)=msqtmp(0)+xmsqjk
            if ((j.eq.0).and.((k.gt.0).or.((k.lt.0).and.
     &           (k.ne.-5)))) then
               bwgttmp(0)=bwgttmp(0)+xmsqjk
            endif
            if ((k.eq.0).and.((j.gt.0).or.((j.lt.0).and.
     &           (j.ne.-5)))) then
               bwgttmp(0)=bwgttmp(0)+xmsqjk
            endif
         endif
cz // end fill index 0

         if (currentPDF .eq. 0) then
           xmsq_bypart(0,sgnj,sgnk)=xmsq_bypart(0,sgnj,sgnk)+xmsqjk

c--->>>P.S. save weight
            if(creategrid.and.bin)then
               currentWeight = msq(j,k)        !*fx1(j)*fx2(k)
               call collectRWeights(0,j,k,currentWeight, realContrib)
               call collectRReference(0,j,k,xmsqjk, procReference)
               if(debug222)then
                  print*,j,k,currentWeight
                  print*,"x1= ",xx1," x2= ",xx2," Q^2= ",facscale**2
                  print*," dipole = ",0
               endif
            endif
c--->>>P.S. end

         endif
         do nd=1,ndmax
	   if (dynamicscale) then	   
             xmsqjk=dipfx1(nd,j)*dipfx2(nd,k)*(-msqc(nd,j,k))
	   else
             xmsqjk=fx1(j)*fx2(k)*(-msqc(nd,j,k))
	   endif
           xmsq(nd)=xmsq(nd)+xmsqjk
           if (currentPDF .eq. 0) then
             xmsq_bypart(nd,sgnj,sgnk)=xmsq_bypart(nd,sgnj,sgnk)+xmsqjk
c--->>>P.S. save weight
             if(creategrid.and.bin)then
                 currentWeight = (-msqc(nd,j,k))  !*fx1(j)*fx2(k)
                 call collectRWeights(nd,j,k,currentWeight,realContrib)
                 call collectRReference(nd,j,k,xmsqjk, procReference)
                 if(debug222)then
                    print*,j,k,currentWeight
                    print*,"x1= ",xx1," x2= ",xx2," Q^2= ",facscale**2
                    print*," dipole = ",nd
                 endif
             endif
c--->>>P.S. end
           endif
         enddo
         
      endif
 20   continue

      enddo
      enddo

c--- loop over all PDF error sets, if necessary
      if (PDFerrors) then
        do nd=0,ndmax 
          PDFxsec_nd(currentPDF,nd)=xmsq(nd)
        enddo
        currentPDF=currentPDF+1
        if (currentPDF .le. maxPDFsets) goto 777
c--- reset xmsq to the central PDF values
        do nd=0,ndmax 
          xmsq(nd)=PDFxsec_nd(0,nd)
        enddo
      endif    

      realint=0d0
      xint=0d0
cz
      bwgt=0d0
cz //

      valsum=0d0 ! running total of weights at this point

c--- zero out temporary histograms
      if (bin) call zerorealhistos

c---trial with weight of real alone
c---first set up all dipole contributions
c---this is the value of integral including subtractions
      do nd=0,ndmax
        xmsq(nd)=xmsq(nd)*flux*pswt/BrnRat
c--->>>P.S. multiply by totalFactor
        if(creategrid.and.bin)then      ! P.S. scale with factor
           totalFactor = flux*pswt/BrnRat
           call scaleRContribution(nd,totalFactor,realContrib)
           call scaleRContribution(nd,totalFactor,procReference)
        endif
c--->>>P.S.end
        if (creatent) then
          wt_gg=xmsq_bypart(nd,0,0)*wgt*flux*pswt/BrnRat/dfloat(itmx)
          wt_gq=(xmsq_bypart(nd,+1,0)+xmsq_bypart(nd,-1,0)
     .          +xmsq_bypart(nd,0,+1)+xmsq_bypart(nd,0,-1)
     .          )*wgt*flux*pswt/BrnRat/dfloat(itmx)
          wt_qq=(xmsq_bypart(nd,+1,+1)+xmsq_bypart(nd,-1,-1)
     .          )*wgt*flux*pswt/BrnRat/dfloat(itmx)
          wt_qqb=(xmsq_bypart(nd,+1,-1)+xmsq_bypart(nd,-1,+1)
     .          )*wgt*flux*pswt/BrnRat/dfloat(itmx)
        endif
        failed=.false.
        
c--- if this dipole has no contribution, go to end of loop
c        if (xmsq(nd) .eq. 0d0) goto 997         
         
        if (nd .eq. 0) then
c---if there's no real contribution, record the event as failing to pass cuts
          if (xmsq(nd) .eq. 0d0) then
             failed=.true.
             goto 996
          endif
        else
c--- if this dipole has no contribution, go to end of loop
          if (xmsq(nd) .eq. 0d0) goto 997         
c---check whether each counter-event passes the cuts
          do j=1,mxpart
          do k=1,4
          q(j,k)=ptilde(nd,j,k)
          enddo
          enddo

c          write(6,*)
c          write(6,*) 'Dipole nd=',nd
          if (incldip(nd)) incldip(nd)=includedipole(nd,q)
          if (incldip(nd) .eqv. .false.) failed=.true.
        endif

 996    if (failed) then
          if (nd .eq. 0) then
            ncutzero=ncutzero+1
            ntotzero=ntotzero+1
          endif
          call dotem(nvec,p,s)
          xmsq(nd)=0d0
          goto 997         
        endif

c---if it does, add to total
        xint=xint+xmsq(nd)
        do j=-1,1
        do k=-1,1
          lord_bypart(j,k)=lord_bypart(j,k)+
     .         wgt*flux*pswt*xmsq_bypart(nd,j,k)/BrnRat
        enddo
        enddo

        val=xmsq(nd)*wgt
        val2=xmsq(nd)**2*wgt
	
	valsum=valsum+val

c     P.S. multiply by totalFactor
        if(creategrid.and.bin)then      ! P.S. scale with factor
           totalFactor = wgt
           call scaleRContribution(nd,totalFactor,realContrib)
           call scaleRContribution(nd,totalFactor,procReference)
        endif
c     P.S.end

cz Fill bwgt if needed
        if(dabs(msqtmp(nd)).gt.0d0) bwgt=bwgttmp(nd)/msqtmp(nd)
cz //
        
c--- update PDF errors
        if (PDFerrors) then
          do currentPDF=0,maxPDFsets        
          PDFwgt(currentPDF)=
     .       flux*pswt*PDFxsec_nd(currentPDF,nd)/BrnRat*wgt/itmx
          PDFxsec(currentPDF)=PDFxsec(currentPDF)
     .       +PDFwgt(currentPDF)
          enddo           
        endif
                
c---if we're binning, add to histo too
        if (bin) then
          call getptildejet(nd,pjet)
          call dotem(nvec,pjet,s)

c     P.S. multiply by totalFactor
          if(creategrid.and.bin)then   ! P.S. scale with factor
c-------      M.S. for some reason he no longer scales by itmx for val
c-------      so remove this here also.
c             totalFactor = 1d0/dfloat(itmx)
              totalFactor = 1d0
              call scaleRContribution(nd,totalFactor,realContrib)
              call scaleRContribution(nd,totalFactor,procReference)
          endif
c     P.S.end

          if (nd .eq. 0) then
            call nplotter(pjet,val,val2,0)
          else
            call nplotter(pjet,val,val2,1)
          endif

c---- >>>>  P.S.  filling grid with current weight
          if (creategrid.and.bin) then

              do iProcess=0,nSubProcess-1
C                 if ( realContrib(nd,iProcess).ne.0 ) then
C                    write (6,*) "contrib ", realContrib(nd,iProcess)
C                 endif   
                 fillContrib(iProcess) = realContrib(nd,iProcess)
              enddo
c------------  QCD coupling removal -------------------------
              totalFactor = eight*pisq/gsq
              call scaleContribution(totalFactor,fillContrib)
c------------  QCD coupling removal -------------------------
c              wrap = yraptwo(3,4,pjet) ! current observable
              call fill_grid(xx1,xx2,facscale,pjet,fillContrib,1)
              call fill_grid_reference(pjet, val, val2)
              if(debug222)then
                 print*,"dipole = ",nd
                 print*,"x1 = ",xx1," x2 = ",xx2," sc = ",facscale
                 print*," gives cont = ",fillContrib
              endif

              do iProcess=0,nSubProcess-1
                 fillContrib(iProcess) = procReference(nd,iProcess)
              enddo
              call fill_process_reference(pjet,fillContrib)
              if(debug222)then
                 print*,"Reference contribution per process"
                 print*," procRef = ",fillContrib
              endif

          endif 
c---- >>>> P.S. end filling grid

        endif
c---otherwise, skip contribution
 997    continue
      enddo

c 998  continue

c--- add temporary histograms to cumulative totals
      if (bin) call addrealhistos(wgt)

c--- update the maximum weight so far, if necessary
        if (dabs(valsum) .gt. wtmax) then
          wtmax=dabs(valsum)
        endif

      if (realwt) then
      realint=xmsq(0)
      else
      realint=xint
      endif
      xreal=xreal+xint*wgt/dfloat(itmx)
      xreal2=xreal2+(xint*wgt)**2/dfloat(itmx)
      

      return

 999  realint=0d0
      ntotzero=ntotzero+1
 
      return
      
      end















